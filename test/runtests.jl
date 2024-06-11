using Test
using PulseAlgorithm

@testset "All Tests" begin

    include("dijkstra.jl")
    include("s_alpha_pulse.jl")
    include("data_loader.jl")

end