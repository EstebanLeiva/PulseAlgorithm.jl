using Test
using PulseAlgorithm

@testset "All Tests" begin

    include("test_dijkstra.jl")
    include("test_pa_algorithm.jl")
    
end