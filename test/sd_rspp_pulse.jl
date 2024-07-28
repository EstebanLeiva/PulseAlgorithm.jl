using DataStructures
using Distributions

@testset "sd_rspp_pulse Test" begin
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

    α = 0.9
    covariance_dict = DefaultDict{Tuple{Int, Int, Int, Int}, Float64}(0.0) 
    covariance_dict[(6, 4, 4, 5)] = 1.0
    covariance_dict[(6, 1, 1, 7)] = 1.0

    pulse = PA.initialize(G, α, covariance_dict, "s", "e")
    PA.preprocess!(pulse)

    @test pulse.mean_costs == [2.0, 9.0, 1.0, 5.0, 2.0, 2.0, 0.0]
    @test pulse.variance_costs == [0.5, 1.0, 0.5, 5.0, 2.0, 1.0, 0.0]

    optimal_path, time, _ = PA.run_pulse(pulse)
    mean, variance, covariance = PA.get_path_distribution(G, optimal_path, covariance_dict)
    quant = quantile(Normal(mean, √(variance + covariance)), α)

    @test optimal_path == [6, 3, 7]
    @test quant - 3.282 <= 1e-3

end