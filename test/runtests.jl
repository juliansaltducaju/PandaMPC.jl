import Pkg; Pkg.add("Test")
using Test, PandaMPC

using ECOS
using Convex
using Plots

@testset "SimulationB" begin
    include("SimulationBtest.jl")
end
