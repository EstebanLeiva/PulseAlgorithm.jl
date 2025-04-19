using Test
using PulseAlgorithm
const PA = PulseAlgorithm

@testset "All Tests" begin

    include("test_data_loader.jl")
    include("test_dijkstra.jl")
    include("test_pa_algorithm.jl")
    
end