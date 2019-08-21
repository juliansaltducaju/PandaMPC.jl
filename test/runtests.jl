using Test, JuliansMPCforJulia

using ECOS
using Convex
using Plots

@testset "SimulationB" begin
    include("SimulationBtest.jl")
end
