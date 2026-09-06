using SARDatasets
using Test

include("fixture.jl")

@testset verbose = true "SARDatasets.jl" begin
    @time @testset "NISAR reader" begin include("nisar.jl") end
    @time @testset "access layer" begin include("remote.jl") end
    # The bridge needs ImagePairGeometry, which is a weak dependency: it is in `test/Project.toml` where
    # it can be resolved, and skipped where it cannot.
    if Base.find_package("ImagePairGeometry") !== nothing
        @time @testset "ImagePairGeometry bridge" begin include("bridge.jl") end
    else
        @info "skipping the bridge test; ImagePairGeometry is not installed"
    end
    # Transfers a few megabytes from a DAAC and needs Earthdata credentials in `~/.netrc`.
    if get(ENV, "SAR_LIVE_TEST", "") == "1"
        @time @testset "live NISAR granule" begin include("live_nisar.jl") end
    else
        @info "skipping the live granule test; set SAR_LIVE_TEST=1 to run it"
    end
end
