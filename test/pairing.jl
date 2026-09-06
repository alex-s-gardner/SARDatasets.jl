# Relating two acquisitions.
#
# Turning a pair into geometry types is ImagePairGeometry's extension and is tested there. What belongs
# here is the interval itself, which is a property of the two products.

using SARDatasets
using Test

@testset "an acquisition against itself has no interval" begin
    mktempdir() do dir
        s = open_sar(write_fixture_product(joinpath(dir, "a.h5")))
        @test repeat_interval(s, s) == 0.0
    end
end

@testset "the interval spans both epochs" begin
    # Two acquisitions weeks apart carry different epochs, so the interval is the epoch difference plus
    # the azimuth difference. Differencing only the azimuth times would give a number under a day for a
    # 48-day pair — which is what the granule this fixture comes from actually is.
    mktempdir() do dir
        a = open_sar(write_fixture_product(joinpath(dir, "a.h5")))
        later = override(FIXTURE,
                         (:geometry => :epoch) => "seconds since 2025-12-15T00:00:00",
                         (:orbit => :epoch) => "seconds since 2025-12-15T00:00:00")
        b = open_sar(write_fixture_product(joinpath(dir, "b.h5"), later))

        @test repeat_interval(a, b) ≈ 48 * 86400 atol = 1.0
        # Antisymmetric: reversing the pair negates it, which is how a caller detects a reversed pair.
        @test repeat_interval(b, a) ≈ -repeat_interval(a, b)
    end
end

@testset "epoch_offset is the epoch's seconds past midnight" begin
    mktempdir() do dir
        # NISAR's epoch is midnight of the acquisition day, so the offset is zero — the case that makes
        # assuming zero look safe.
        s = open_sar(write_fixture_product(joinpath(dir, "a.h5")))
        @test epoch_offset(s) == 0.0
        @test epoch_offset(s.geometry) == epoch_offset(s)

        # A product on any other epoch, which is why it is computed.
        noon = override(FIXTURE,
                        (:geometry => :epoch) => "seconds since 2025-10-28T12:00:00",
                        (:orbit => :epoch) => "seconds since 2025-10-28T12:00:00")
        @test epoch_offset(open_sar(write_fixture_product(joinpath(dir, "b.h5"), noon))) == 43200.0
    end
end
