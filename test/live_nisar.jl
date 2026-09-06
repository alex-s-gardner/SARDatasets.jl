# Reading a real granule over the network.
#
# Needs Earthdata credentials in `~/.netrc` and transfers a few megabytes, so it runs only when
# `SAR_LIVE_TEST=1`. What it proves that the fixture tests cannot: that a multi-gigabyte product on a
# DAAC opens through the redirect chain, and that only the prefetch window is transferred.
#
# The values asserted are the committed fixture's, which was harvested from this granule with `h5py`.
# So a pass here means the Julia reader and the reference agree on a real product, bitwise.

using SARDatasets
using Test

const LIVE_URL = get(ENV, "SAR_LIVE_URL",
    "https://nisar.asf.earthdatacloud.nasa.gov/NISAR/NISAR_L1_RSLC_BETA_V1/" *
    "NISAR_L1_PR_RSLC_003_170_D_053_7700_SHNA_A_20251028T235201_20251028T235238_X05009_N_P_J_001/" *
    "NISAR_L1_PR_RSLC_003_170_D_053_7700_SHNA_A_20251028T235201_20251028T235238_X05009_N_P_J_001.h5")

mktempdir() do dir
    src = RemoteHTTP(LIVE_URL; dir)
    s = open_sar(src)
    o = orbit(s)
    path = only(readdir(dir; join = true))

    @testset "only the prefetch window was transferred" begin
        apparent = filesize(path)
        # Blocks the filesystem actually allocated: a sparse file's holes occupy none, so this is what
        # crossed the network.
        allocated = 512 * parse(Int, first(split(read(`stat -f "%b" $path`, String))))
        @test apparent > 1_000_000_000            # a real granule, not a stub
        @test allocated <= src.prefetch
        @test allocated < apparent ÷ 100          # far less than one percent of it
        @info "live granule" apparent allocated fraction = allocated / apparent
    end

    @testset "the live read matches the committed fixture bitwise" begin
        g = s.geometry
        @test g.starting_range === gx(FIXTURE.geometry.starting_range)
        @test g.far_range === gx(FIXTURE.geometry.far_range)
        @test g.range_pixel_spacing === gx(FIXTURE.geometry.range_pixel_spacing)
        @test g.wavelength === gx(FIXTURE.geometry.wavelength)
        @test g.prf === gx(FIXTURE.geometry.prf)
        @test g.sensing_start === gx(FIXTURE.geometry.sensing_start)
        @test g.sensing_stop === gx(FIXTURE.geometry.sensing_stop)
        @test g.nlines == FIXTURE.geometry.nlines
        @test g.nsamples == FIXTURE.geometry.nsamples

        @test length(o.time) == FIXTURE.orbit.n
        @test all(o.time[i] === gx(FIXTURE.orbit.time[i]) for i in eachindex(o.time))
        @test all(o.position[i][c] === gx(FIXTURE.orbit.position[i][c])
                  for i in eachindex(o.position), c in 1:3)
        @test all(o.velocity[i][c] === gx(FIXTURE.orbit.velocity[i][c])
                  for i in eachindex(o.velocity), c in 1:3)

        @test s.identification.mission == FIXTURE.identification.mission
        @test s.identification.absolute_orbit == FIXTURE.identification.absolute_orbit
        @test s.identification.start_time == FIXTURE.identification.start_time
    end

    @testset "a prefetch too small to hold the metadata names the knob" begin
        err = try
            open_sar(RemoteHTTP(LIVE_URL; prefetch = 64 * 1024, dir = mktempdir()))
            nothing
        catch e
            sprint(showerror, e)
        end
        @test err !== nothing
        @test occursin("prefetch", err)
        @test occursin("65536", err)
    end
end
