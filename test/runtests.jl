using SARDatasets
using Test

include("fixture.jl")

# Whether the geometry bridge can be tested, which is whether ImagePairGeometry can be loaded at all.
#
# Asked by trying rather than by looking: `Base.find_package` answers about the active project, and
# `Pkg.test` runs in a sandbox that sees the test environment while the package environment is where an
# unregistered weak dependency has to be added — Julia 1.10 refuses to merge the two when a git-sourced
# package appears in both. So the reliable question is whether `import` succeeds here.
const HAS_GEOMETRY = try
    @eval import ImagePairGeometry
    true
catch
    false
end

@testset verbose = true "SARDatasets.jl" begin
    @time @testset "NISAR reader" begin include("nisar.jl") end
    @time @testset "access layer" begin include("remote.jl") end
    if HAS_GEOMETRY
        @time @testset "ImagePairGeometry bridge" begin include("bridge.jl") end
    else
        @info "skipping the bridge test; ImagePairGeometry could not be loaded"
    end
    # Transfers a few megabytes from a DAAC and needs Earthdata credentials in `~/.netrc`.
    if get(ENV, "SAR_LIVE_TEST", "") == "1"
        @time @testset "live NISAR granule" begin include("live_nisar.jl") end
    else
        @info "skipping the live granule test; set SAR_LIVE_TEST=1 to run it"
    end
end
