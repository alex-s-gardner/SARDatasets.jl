# The NISAR path, on a fixture product built the way the test suite builds one, plus a count of how
# many times one `open_slc` opens the file.

using BenchmarkTools, Printf, HDF5, JSON3
using SLCDatasets
using SLCDatasets: SPEED_OF_LIGHT, NisarBackend, read_identification, read_geometry, read_orbit,
                   nisar_band, nisar_product_type, default_frequency, orbit

const WORK = get(ENV, "SLC_BENCH_DIR", "/tmp/slcbench")
mkpath(WORK)
const H5 = joinpath(WORK, "nisar_fixture.h5")

# A product at a real RSLC's grid size, with the orbit vector count a real one carries.
function write_nisar(path; nlines = 27_000, nsamples = 12_000, nsv = 1_500)
    isfile(path) && return path
    h5open(path, "w") do h
        p, ip = "science/LSAR/RSLC", "science/LSAR/identification"
        h["$ip/missionId"] = "NISAR"
        h["$ip/productType"] = "RSLC"
        h["$ip/absoluteOrbitNumber"] = Int32(1234)
        h["$ip/orbitPassDirection"] = "Descending"
        h["$ip/lookDirection"] = "Left"
        h["$ip/zeroDopplerStartTime"] = "2025-10-28T12:00:00.123456789"
        h["$ip/zeroDopplerEndTime"] = "2025-10-28T12:00:28.987654321"
        h["$ip/boundingPolygon"] = "POLYGON EMPTY"
        h["$ip/listOfFrequencies"] = ["A"]
        h["$p/swaths/zeroDopplerTime"] = collect(range(43200.0, 43228.0; length = nlines))
        h["$p/swaths/zeroDopplerTimeSpacing"] = 28.0 / (nlines - 1)
        h["$p/swaths/frequencyA/slantRange"] = collect(range(8.0e5, 9.0e5; length = nsamples))
        h["$p/swaths/frequencyA/slantRangeSpacing"] = 1.0e5 / (nsamples - 1)
        h["$p/swaths/frequencyA/processedCenterFrequency"] = 1.2575e9
        h["$p/metadata/orbit/time"] = collect(range(43100.0, 43300.0; length = nsv))
        h["$p/metadata/orbit/position"] = randn(3, nsv) .* 1e6
        h["$p/metadata/orbit/velocity"] = randn(3, nsv) .* 1e3
        h["$p/metadata/orbit/interpMethod"] = "Hermite"
        h["$p/metadata/orbit/orbitType"] = "POE"
        write_attribute(h["$p/swaths/zeroDopplerTime"], "units",
                        "seconds since 2025-10-28 00:00:00")
        write_attribute(h["$p/metadata/orbit/time"], "units",
                        "seconds since 2025-10-28 00:00:00")
    end
    return path
end

write_nisar(H5)
@printf("fixture: %.1f MB\n\n", filesize(H5) / 1e6)

b = NisarBackend(H5, "LSAR", "RSLC", "A")

function show_case(name, f)
    f()
    t = @benchmark $f() samples = 20 evals = 1 seconds = 25
    @printf("%-36s %10s %10d %10s\n", name, BenchmarkTools.prettytime(minimum(t.times)),
            minimum(t.allocs), BenchmarkTools.prettymemory(minimum(t.memory)))
end

@printf("%-36s %10s %10s %10s\n", "case", "time", "allocs", "memory")
println("-"^70)
show_case("open_slc (NISAR)", () -> open_slc(H5))
show_case("open_slc + orbit", () -> orbit(open_slc(H5)))
println()
show_case("nisar_band", () -> nisar_band(H5))
show_case("nisar_product_type", () -> nisar_product_type(H5, "LSAR"))
show_case("default_frequency", () -> default_frequency(H5, "LSAR"))
show_case("read_identification", () -> read_identification(b))
show_case("read_geometry", () -> read_geometry(b))
show_case("read_orbit", () -> read_orbit(b))
println()
show_case("bare h5open/close", () -> h5open(identity, H5, "r"))

# `read_geometry` reads the whole 27k-element azimuth axis and 12k-element range axis to take their
# first and last elements. Reading just the endpoints is the alternative.
println()
show_case("read whole zeroDopplerTime (27k)",
          () -> h5open(h -> read(h["science/LSAR/RSLC/swaths/zeroDopplerTime"]), H5, "r"))
show_case("read its two endpoints only",
          () -> h5open(H5, "r") do h
              d = h["science/LSAR/RSLC/swaths/zeroDopplerTime"]
              n = length(d)
              (d[1], d[n], n)
          end)
