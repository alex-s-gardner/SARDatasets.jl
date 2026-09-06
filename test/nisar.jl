# The NISAR reader against the committed fixture.
#
# Floats are compared with `===` against the fixture's hex literals, so these assert bit-exact
# agreement with what `h5py` read from the granule, not agreement to a printed precision.

using SARDatasets
using SARDatasets: LookSide, LookLeft, LookRight, SPEED_OF_LIGHT, parse_cf_epoch, NisarBackend, nisar_band,
           nisar_product_type, GEOCODED_TYPES
using Dates
using Test

const FX = FIXTURE

mktempdir() do dir
    path = write_fixture_product(joinpath(dir, "fixture_rslc.h5"))
    s = open_sar(path)

    @testset "identification" begin
        id = s.identification
        @test id.mission == FX.identification.mission
        @test id.product_type == FX.identification.product_type
        @test id.absolute_orbit == FX.identification.absolute_orbit
        @test id.pass_direction == FX.identification.pass_direction
        @test id.look_direction == FX.identification.look_direction
        # Kept as strings: the product records nanoseconds, which a `DateTime` would silently drop.
        @test id.start_time == FX.identification.start_time
        @test id.stop_time == FX.identification.stop_time
        @test start_datetime(id) == DateTime("2025-10-28T23:52:01.000")
        @test stop_datetime(id) == DateTime("2025-10-28T23:52:38.999")
    end

    @testset "geometry is bitwise" begin
        g = s.geometry
        @test g.starting_range === gx(FX.geometry.starting_range)
        @test g.far_range === gx(FX.geometry.far_range)
        @test g.range_pixel_spacing === gx(FX.geometry.range_pixel_spacing)
        @test g.wavelength === gx(FX.geometry.wavelength)
        @test g.prf === gx(FX.geometry.prf)
        @test g.sensing_start === gx(FX.geometry.sensing_start)
        @test g.sensing_stop === gx(FX.geometry.sensing_stop)
        @test g.nlines == FX.geometry.nlines
        @test g.nsamples == FX.geometry.nsamples
        @test g.look_side == (lowercase(FX.identification.look_direction) == "left" ?
                              LookLeft : LookRight)
        @test nlines(s) == g.nlines
        @test nsamples(s) == g.nsamples
    end

    @testset "wavelength is c over the center frequency" begin
        # The reference derives the wavelength this way, so the division must not be reordered.
        fc = SPEED_OF_LIGHT / gx(FX.geometry.wavelength)
        @test s.geometry.wavelength === SPEED_OF_LIGHT / fc
    end

    @testset "orbit is bitwise" begin
        o = orbit(s)
        @test length(o.time) == FX.orbit.n
        @test o.interp_method == FX.orbit.interp_method
        @test o.kind == FX.orbit.kind
        @test all(o.time[i] === gx(FX.orbit.time[i]) for i in eachindex(o.time))
        # The file stores (N, 3); a transposed read would swap x and z here and pass a norm check.
        @test all(o.position[i][c] === gx(FX.orbit.position[i][c])
                  for i in eachindex(o.position), c in 1:3)
        @test all(o.velocity[i][c] === gx(FX.orbit.velocity[i][c])
                  for i in eachindex(o.velocity), c in 1:3)
    end

    @testset "the orbit is read once and held" begin
        s2 = open_sar(path)
        @test s2.orbit === nothing
        first_read = orbit(s2)
        @test s2.orbit !== nothing
        @test orbit(s2) === first_read
    end

    @testset "one epoch for both clocks" begin
        # The azimuth times and the state vector times share an epoch in this format. A consumer
        # converting between them relies on it, so a product where they diverge must be visible.
        @test s.geometry.epoch == orbit(s).epoch
        @test s.geometry.epoch == parse_cf_epoch(FX.geometry.epoch)
    end

    @testset "the orbit brackets the acquisition" begin
        o = orbit(s)
        @test first(o.time) <= s.geometry.sensing_start
        @test s.geometry.sensing_stop <= last(o.time)
    end

    @testset "backend paths" begin
        b = s.backend
        @test b isa NisarBackend
        @test b.band == FX.band
        @test b.product_type == FX.product_type
        @test b.frequency == FX.frequency
        @test nisar_band(path) == FX.band
        @test nisar_product_type(path, FX.band) == FX.product_type
    end

    @testset "show" begin
        @test occursin("RSLC", sprint(show, s))
        long = sprint(show, MIME"text/plain"(), s)
        @test occursin("left", long)
        @test occursin("35 state vectors", sprint(show, MIME"text/plain"(), s))
    end
end

@testset "parse_cf_epoch" begin
    @test parse_cf_epoch("seconds since 2025-10-28T00:00:00") == DateTime(2025, 10, 28)
    @test parse_cf_epoch("seconds since 2025-10-28 00:00:00") == DateTime(2025, 10, 28)
    @test parse_cf_epoch("seconds since 2025-10-28T00:00:00Z") == DateTime(2025, 10, 28)
    # Rescaling a different unit would corrupt every time by a constant factor, so it throws.
    @test_throws "only seconds are supported" parse_cf_epoch("days since 2025-10-28T00:00:00")
    @test_throws "cannot read a reference epoch" parse_cf_epoch("2025-10-28")
end

@testset "unreadable inputs throw" begin
    @test_throws "is not a readable file" open_sar(joinpath(@__DIR__, "no_such_file.h5"))
    mktempdir() do dir
        plain = joinpath(dir, "plain.txt")
        write(plain, "not hdf5")
        @test_throws "is not an HDF5 file" open_sar(plain)

        empty = joinpath(dir, "empty.h5")
        h5open(empty, "w") do h
            h["unrelated"] = 1
        end
        @test_throws "not a NISAR-format product" open_sar(empty)
    end
end

@testset "a geocoded product is rejected, not misread" begin
    # A GSLC carries map-projected grids and no slant-range axis. Reading one as if it had a swath
    # would fail on a missing group; the reader says why instead.
    mktempdir() do dir
        path = joinpath(dir, "gslc.h5")
        h5open(path, "w") do h
            h["science/LSAR/identification/productType"] = "GSLC"
            h["science/LSAR/identification/absoluteOrbitNumber"] = Int32(1)
            h["science/LSAR/identification/listOfFrequencies"] = ["A"]
            h["science/LSAR/GSLC/grids/frequencyA/x"] = [0.0, 1.0]
        end
        @test "GSLC" in GEOCODED_TYPES
        @test_throws "carries no slant-range/azimuth geometry" open_sar(path)
    end
end
