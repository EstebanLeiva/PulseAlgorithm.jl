using PulseAlgorithm
using Test
using DataStructures
using Distributions

@testset "s_alpha_pulse" begin
    G = Graph(Dict{Int, Node}(), Dict{String, Int}())

    create_node!(G, "1")
    create_node!(G, "2")
    create_node!(G, "3")
    create_node!(G, "4")
    create_node!(G, "5")
    create_node!(G, "s") #6
    create_node!(G, "e") #7

    add_link!(G, "s", "1", 2.0, 2.0, 3.0)
    add_link!(G, "1", "e", 3.0, 2.0, 0.5)
    add_link!(G, "s", "2", 3.0, 2.0, 1.0)
    add_link!(G, "2", "e", 5.0, 9.0, 1.0)
    add_link!(G, "s", "3", 2.0, 1.0, 0.5)
    add_link!(G, "3", "e", 4.0, 1.0, 0.5)
    add_link!(G, "s", "4", 1.0, 2.0, 3.0)
    add_link!(G, "4", "5", 1.0, 3.0, 3.0)
    add_link!(G, "5", "e", 1.0, 2.0, 2.0)

    @info "Graph created"

    T = 10.0
    α = 0.9
    covariance_dict = DefaultDict{Tuple{Int, Int, Int, Int}, Float64}(0.0) 
    covariance_dict[(6, 4, 4, 5)] = 1.0
    covariance_dict[(6, 1, 1, 7)] = 1.0

    pulse = create_SPulseGraph(G, α, covariance_dict, "s", "e", T)
    preprocess!(pulse)

    @test pulse.minimum_costs == [3.0, 5.0, 4.0, 2.0, 1.0, 3.0, 0.0]
    @test pulse.mean_costs == [2.0, 9.0, 1.0, 5.0, 2.0, 2.0, 0.0]
    @test pulse.variance_costs == [0.5, 1.0, 0.5, 5.0, 2.0, 1.0, 0.0]

    optimal_path, cost, _ = run_pulse(pulse)
    mean, variance, covariance = get_path_distribution(G, optimal_path, covariance_dict)
    reliability = cdf(Normal(mean, √(variance + covariance)), T)

    @test optimal_path == [6, 1, 7]
    @test cost == 5.0
    @test reliability >= 0.99

end