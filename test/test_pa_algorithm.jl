using Distributions
using PulseAlgorithm: Graph, Parameters, Problem, Pulse, preprocess!, run_pulse!, DefaultDict

@testset "S-aRP Test" begin

    function get_path_distribution(graph::Graph, path::Vector{Int}, cov_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64})
        mean = 0.0
        variance = 0.0
        covariance_term = 0.0
        for i in 1:length(path)-1
            mean += graph.nodes[path[i]].links[path[i+1]].random["time"]["mean"]
            variance += graph.nodes[path[i]].links[path[i+1]].random["time"]["variance"]
            for ii in i + 1:length(path)-1
                covariance_term += 2*cov_dict[(path[i], path[i+1], path[ii], path[ii+1])]
            end
        end
        return mean, variance, covariance_term
    end

    function info_update(graph::Graph,
                         current_node::Int, 
                         reachable_node::Int, 
                         path::Vector{Int},
                         deterministic_info::Dict{String, Float64}, 
                         random_info::Dict{String, Dict{String, Float64}})
        deterministic_info = copy(deterministic_info) #TODO: check if this copying can be done automatically in some way
        random_info = copy(random_info)

        deterministic_info["cost"] += graph.nodes[current_node].links[reachable_node].deterministic["cost"]
        random_info["time"]["mean"] += graph.nodes[current_node].links[reachable_node].random["time"]["mean"]
        random_info["time"]["variance"] += graph.nodes[current_node].links[reachable_node].random["time"]["variance"]
        
        n = length(path)
        if n > 1
            last_node = path[end]
            current_node = reachable_node
            covariance_sum = 0.0
            for i in 1:n-1
                covariance_sum += 2 * graph.covariance["time"][(path[i], path[i+1], last_node, current_node)]
            end
        else
            covariance_sum =  0.0
        end

        random_info["time"]["covariance"] += covariance_sum

        return deterministic_info, random_info
    end

    function prune_feasibility(pulse_alg::Pulse, 
                               current_node::Int, 
                               current_path::Vector{Int},
                               deterministic_info::Dict{String, Float64},
                               random_info::Dict{String, Dict{String, Float64}})
        println(current_path)
        pass = false
        mean = random_info["time"]["mean"] + pulse_alg.prep_random_costs["time"]["mean"][current_node]
        variance = random_info["time"]["variance"] + pulse_alg.prep_random_costs["time"]["variance"][current_node]
        dist = Normal(mean, √variance)
        prob = cdf(dist, pulse_alg.problem.constants["T_max"])
        if pulse_alg.problem.constants["T_max"] >= mean && prob < pulse_alg.problem.constants["alpha"]
            pass = true
        elseif pulse_alg.problem.constants["T_max"] < mean && pulse_alg.problem.constants["alpha"] > 0.5
            pass = true
        end
        return pass
    end

    function prune_bounds(pulse_alg::Pulse, 
                          current_node::Int, 
                          current_path::Vector{Int},
                          deterministic_info::Dict{String, Float64},
                          random_info::Dict{String, Dict{String, Float64}})
        pass = true
        if deterministic_info["cost"] + pulse_alg.prep_deterministic_costs["cost"][current_node] <= pulse_alg.current_optimal_objective
            if current_node == pulse_alg.problem.target_node
                pulse_alg.current_optimal_objective = deterministic_info["cost"]
                new_path = copy(current_path)
                push!(new_path, current_node)
                pulse_alg.current_optimal_path = new_path
            end
            pass = false
        end
        return pass
    end

    function exploration_order(pulse_alg::Pulse, node::Int)
        return pulse_alg.prep_deterministic_costs["cost"][node]
    end
    
    function pulse_score(pulse_alg::Pulse, 
                         current_path::Vector{Int},
                         deterministic_info::Dict{String, Float64},
                         random_info::Dict{String, Dict{String, Float64}})
        return pulse_alg.prep_deterministic_costs["cost"][current_path[end]]
    end
        
    pruning_functions = [prune_bounds, prune_feasibility]

    cov = Dict("time" => DefaultDict{Tuple{Int, Int, Int, Int}, Float64}(0.0))
    G = Graph(Dict{Int, Node}(), Dict{String, Int}(), cov)

    create_node!(G, "1")
    create_node!(G, "2")
    create_node!(G, "3")
    create_node!(G, "4")
    create_node!(G, "5")
    create_node!(G, "s") #6
    create_node!(G, "e") #7

    add_link!(G, "s", "1", Dict("cost" => 2.0), Dict("time" => Dict("mean" => 2.0, "variance" => 3.0)))
    add_link!(G, "1", "e", Dict("cost" => 3.0), Dict("time" => Dict("mean" => 2.0, "variance" => 0.5)))
    add_link!(G, "s", "2", Dict("cost" => 3.0), Dict("time" => Dict("mean" => 2.0, "variance" => 1.0)))
    add_link!(G, "2", "e", Dict("cost" => 5.0), Dict("time" => Dict("mean" => 9.0, "variance" => 1.0)))
    add_link!(G, "s", "3", Dict("cost" => 2.0), Dict("time" => Dict("mean" => 1.0, "variance" => 0.5)))
    add_link!(G, "3", "e", Dict("cost" => 4.0), Dict("time" => Dict("mean" => 1.0, "variance" => 0.5)))
    add_link!(G, "s", "4", Dict("cost" => 1.0), Dict("time" => Dict("mean" => 2.0, "variance" => 3.0)))
    add_link!(G, "4", "5", Dict("cost" => 1.0), Dict("time" => Dict("mean" => 3.0, "variance" => 3.0)))
    add_link!(G, "5", "e", Dict("cost" => 1.0), Dict("time" => Dict("mean" => 2.0, "variance" => 2.0)))

    G.covariance["time"][(6, 4, 4, 5)] = 1.0
    G.covariance["time"][(6, 1, 1, 7)] = 1.0
    
    # Parameters
    deterministic_weights = ["cost"]
    random_weights = Dict("time" => ["mean", "variance", "covariance"])
    prep_deterministic_weights = ["cost"]
    prep_random_weights = Dict("time" => ["mean", "variance"])
    params = Parameters(1000000, false, exploration_order, deterministic_weights, random_weights, prep_deterministic_weights, prep_random_weights)

    # Problem
    constants = Dict("T_max" => 10.0, "alpha" => 0.9)
    problem = Problem(G, 6, 7, true, constants)

    pulse = Pulse(problem, params)

    preprocess!(pulse)

    @test pulse.prep_deterministic_costs["cost"] == [3.0, 5.0, 4.0, 2.0, 1.0, 3.0, 0.0]
    @test pulse.prep_random_costs["time"]["mean"] == [2.0, 9.0, 1.0, 5.0, 2.0, 2.0, 0.0]
    @test pulse.prep_random_costs["time"]["variance"] == [0.5, 1.0, 0.5, 5.0, 2.0, 1.0, 0.0]

    run_pulse!(pulse, 
               info_update, 
               pruning_functions, 
               pulse_score)

    optimal_path = pulse.optimal_path
    optimal_cost = pulse.optimal_objective

    mean, variance, covariance = get_path_distribution(G, optimal_path, cov["time"])
    reliability = cdf(Normal(mean, √(variance + covariance)), constants["T_max"])

    @test optimal_path == [6, 1, 7]
    @test optimal_cost == 5.0
    @test reliability >= 0.99

end