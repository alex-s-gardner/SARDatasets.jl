# Rebuilding a NISAR-layout product from the committed fixture.
#
# The fixture records what `h5py` read from a real granule. Writing it back out in the same group
# layout gives the reader a product to open, so the whole path — root-path probe, product-type
# dispatch, group paths, epoch parsing, the (3, N) orbit transpose — is exercised without a multi-
# gigabyte file or a network.

using HDF5
using JSON3

const FIXTURE = JSON3.read(read(joinpath(@__DIR__, "reference", "nisar_metadata.json"), String))

# Every float in the fixture carries a hex literal beside its decimal, and the hex is what the tests
# compare against: a decimal round-trip through JSON is not guaranteed to preserve the last bit.
gx(v) = parse(Float64, v.hex)

"""
    write_fixture_product(path, fx = FIXTURE) -> String

Write `fx` as a NISAR-layout HDF5 product at `path`.

Only the datasets the reader reads are written, with the group names and the `units` attributes a real
product carries.
"""
function write_fixture_product(path::AbstractString, fx = FIXTURE)
    band, product, freq = fx.band, fx.product_type, fx.frequency
    id, geom, orb = fx.identification, fx.geometry, fx.orbit

    n = orb.n
    time = [gx(v) for v in orb.time]
    # The fixture stores state vectors as n rows of 3, which is the file's own (N, 3) layout; HDF5.jl
    # writes a (3, n) Julia array to that shape.
    pos = [gx(orb.position[i][c]) for c in 1:3, i in 1:n]
    vel = [gx(orb.velocity[i][c]) for c in 1:3, i in 1:n]

    h5open(path, "w") do h
        p = "science/$band/$product"
        ip = "science/$band/identification"

        h["$ip/missionId"] = id.mission
        h["$ip/productType"] = id.product_type
        h["$ip/absoluteOrbitNumber"] = Int32(id.absolute_orbit)
        h["$ip/orbitPassDirection"] = id.pass_direction
        h["$ip/lookDirection"] = id.look_direction
        h["$ip/zeroDopplerStartTime"] = id.start_time
        h["$ip/zeroDopplerEndTime"] = id.stop_time
        h["$ip/boundingPolygon"] = "POLYGON EMPTY"
        h["$ip/listOfFrequencies"] = [freq]

        # The azimuth axis is regenerated from its ends and the line count so the file carries a full
        # `zeroDopplerTime` vector of the product's real length, as `nlines` is read from its size.
        t0, t1, nl = gx(geom.sensing_start), gx(geom.sensing_stop), geom.nlines
        h["$p/swaths/zeroDopplerTime"] = collect(range(t0, t1; length = nl))
        h["$p/swaths/zeroDopplerTimeSpacing"] = 1 / gx(geom.prf)

        r0, r1, ns = gx(geom.starting_range), gx(geom.far_range), geom.nsamples
        h["$p/swaths/frequency$freq/slantRange"] = collect(range(r0, r1; length = ns))
        h["$p/swaths/frequency$freq/slantRangeSpacing"] = gx(geom.range_pixel_spacing)
        h["$p/swaths/frequency$freq/processedCenterFrequency"] =
            SPEED_OF_LIGHT / gx(geom.wavelength)
        h["$p/swaths/frequency$freq/listOfPolarizations"] = ["HH"]

        h["$p/metadata/orbit/time"] = time
        h["$p/metadata/orbit/position"] = pos
        h["$p/metadata/orbit/velocity"] = vel
        h["$p/metadata/orbit/interpMethod"] = orb.interp_method
        h["$p/metadata/orbit/orbitType"] = orb.kind

        write_attribute(h["$p/swaths/zeroDopplerTime"], "units", geom.epoch)
        write_attribute(h["$p/metadata/orbit/time"], "units", orb.epoch)
    end
    return path
end
