using SAR
using Test

include("fixture.jl")

@testset verbose = true "SAR.jl" begin
    @time @testset "NISAR reader" begin include("nisar.jl") end
    @time @testset "access layer" begin include("remote.jl") end
    # Transfers a few megabytes from a DAAC and needs Earthdata credentials in `~/.netrc`.
    if get(ENV, "SAR_LIVE_TEST", "") == "1"
        @time @testset "live NISAR granule" begin include("live_nisar.jl") end
    else
        @info "skipping the live granule test; set SAR_LIVE_TEST=1 to run it"
    end
end
