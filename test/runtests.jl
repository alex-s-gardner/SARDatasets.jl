using SAR
using Test

include("fixture.jl")

@testset verbose = true "SAR.jl" begin
    @time @testset "NISAR reader" begin include("nisar.jl") end
end
