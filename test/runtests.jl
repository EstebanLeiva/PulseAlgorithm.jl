using Test
using PulseAlgorithm
const PA = PulseAlgorithm

@testset "All Tests" begin

    include("dijkstra.jl")
    include("s_alpha_pulse.jl")
    include("data_loader.jl")
    include("util.jl")
    include("erspa_star.jl")
    include("sd_rspp_pulse.jl")

end