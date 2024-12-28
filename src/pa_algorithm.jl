struct Parameters
    max_pulse_depth::Int
    path_completion::Bool
    exploration_order::Function
    deterministic_weights::Vector{String}
    random_weights::Dict{String, Vector{String}}
    prep_deterministic_weights::Vector{String}
    prep_random_weights::Dict{String, Vector{String}}
end

function Parameters(json_dir::String)
    #TODO: load parameters from json 
    return 1
end

mutable struct Pulse
    # Problem information
    const problem::Problem
    # Pulse Parameters
    const parameters::Parameters
    # Preprocessing information
    const prep_deterministic_costs::Dict{String, Vector{Float64}}
    const prep_random_costs::Dict{String, Dict{String, Vector{Float64}}}
    # Optimal path information
    optimal_path::Vector{Int}
    optimal_objective::Float64
    # Pruning strategies
    current_optimal_path::Vector{Int}
    current_optimal_objective::Float64
    const dominance::Dict{Tuple{Int, Int}, PriorityQueue{Tuple{Float64, Float64}, Float64}}
    # Acceleration strategies
    const pulse_queue::PriorityQueue{Tuple{Vector{Int}, Dict{String, Float64}, Dict{String, Dict{String, Float64}}}, Float64}
    # Instance information
    const instance_info::Dict{String, Int} #TODO: Add key initialization to 0
end

function Pulse(problem::Problem, parameters::Parameters)
    det_keys, random_keys = get_link_keys(problem.graph)
    if !(parameters.prep_deterministic_weights ⊆ det_keys)
        error("The deterministic weights for preprocessing are not valid. Choose a subset of the following: $det_keys")
    elseif dict_not_subset(parameters.prep_random_weights, random_keys)
        error("The random weights for preprocessing are not valid.")
    end
    return Pulse(problem, 
                 parameters, 
                 Dict{String, Vector{Float64}}(), 
                 Dict{String, Dict{String, Vector{Float64}}}(), 
                 Vector{Int}(), 
                 Inf, 
                 Vector{Int}(), 
                 Inf, 
                 Dict{Tuple{Int, Int}, PriorityQueue{Tuple{Float64, Float64}, Float64}}(), 
                 PriorityQueue{Tuple{Vector{Int}, Dict{String, Float64}, Dict{String, Dict{String, Float64}}}, Float64}(), 
                 Dict{String, Int}())
end

function preprocess!(pulse_alg::Pulse)
    for cost in pulse_alg.parameters.prep_deterministic_weights
        pulse_alg.prep_deterministic_costs[cost] = dijkstra(pulse_alg.problem.graph, pulse_alg.problem.target_node, cost, pulse_alg.parameters.path_completion) 
        if pulse_alg.prep_deterministic_costs[cost][pulse_alg.problem.source_node] == Inf
            error("The source node is not reachable from the target node")
        end
        #TODO: Add path completion
    end
    for random_variable in keys(pulse_alg.parameters.prep_random_weights)
        pulse_alg.prep_random_costs[random_variable] = Dict{String, Vector{Float64}}()
        for cost in pulse_alg.parameters.prep_random_weights[random_variable]
            pulse_alg.prep_random_costs[random_variable][cost] = dijkstra(pulse_alg.problem.graph, pulse_alg.problem.target_node, random_variable, cost, pulse_alg.parameters.path_completion)
            if pulse_alg.prep_random_costs[random_variable][cost][pulse_alg.problem.source_node] == Inf
                error("The source node is not reachable from the target node")
            end
        end
        #TODO: Add path completion
    end 
end

function propagate_pulse!(pulse_alg::Pulse, 
                         current_node::Int,
                         deterministic_info::Dict{String, Float64},
                         random_info::Dict{String, Dict{String, Float64}},
                         current_path::Vector{Int}, 
                         current_depth::Int,
                         pruning_functions::Vector{Function}, 
                         info_update::Function, 
                         pulse_score::Function)
    pass = true
    for pruning_function in pruning_functions
        if pruning_function(pulse_alg, current_node, current_path, deterministic_info, random_info)
            pass = false
            break
        end
    end
    if pass
        push!(current_path, current_node)
        link_dict = pulse_alg.problem.graph.nodes[current_node].links
        if current_path[end] ≠ pulse_alg.problem.target_node
            if current_depth < pulse_alg.parameters.max_pulse_depth
                ordered_reachable_nodes = order_nodes(pulse_alg, link_dict, pulse_alg.parameters.exploration_order)
                for reachable_node in ordered_reachable_nodes
                    if reachable_node ∉ current_path
                        new_path = copy(current_path)
                        new_deterministic_info, new_random_info = info_update(pulse_alg.problem.graph,
                                                                              current_node, 
                                                                              reachable_node, 
                                                                              new_path,
                                                                              deterministic_info, 
                                                                              random_info) #TODO: check if this is a copy of the info
                        propagate_pulse!(pulse_alg, 
                                        reachable_node, 
                                        new_deterministic_info, 
                                        new_random_info, 
                                        new_path, 
                                        current_depth + 1, 
                                        pruning_functions, 
                                        info_update, 
                                        pulse_score)
                    end
                end
            else
                score = pulse_score(pulse_alg, current_path, deterministic_info, random_info)
                enqueue!(pulse_alg.pulse_queue, (current_path, deterministic_info, random_info), score)
            end
        end
    end
end

function run_pulse!(pulse_alg::Pulse,
                   info_update::Function,
                   pruning_functions::Vector{Function},
                   pulse_score::Function,
                   init_optimal_path::Vector{Int} = Vector{Int}(), 
                   init_objective::Float64 = Inf)
    path = Vector{Int}()
    pulse_alg.current_optimal_path = init_optimal_path
    pulse_alg.current_optimal_objective = init_objective
    deterministic_info, random_info = init_info(pulse_alg.parameters.deterministic_weights, pulse_alg.parameters.random_weights)
    propagate_pulse!(pulse_alg, 
                    pulse_alg.problem.source_node, 
                    deterministic_info, 
                    random_info, 
                    path, 
                    0, 
                    pruning_functions, 
                    info_update, 
                    pulse_score)
                    
    while !isempty(pulse_alg.pulse_queue)
        path_to_explore, deterministic_info, random_info = dequeue!(pulse_alg.pulse_queue)
        link_dict = pulse_alg.problem.graph.nodes[path_to_explore[end]].links
        ordered_reachable_nodes = order_nodes(link_dict, pulse_alg.parameters.exploration_order)
        for reachable_node in ordered_reachable_nodes
            if reachable_node ∉ path_to_explore
                inside_path = copy(path_to_explore)
                new_deterministic_info, new_random_info = info_update(path_to_explore[end], 
                                                                      reachable_node, 
                                                                      deterministic_info, 
                                                                      random_info) #TODO: check if this is a copy of the info
                propagate_pulse!(pulse_alg, 
                                reachable_node, 
                                new_deterministic_info, 
                                new_random_info, 
                                inside_path, 
                                length(inside_path), 
                                pruning_functions, 
                                info_update, 
                                pulse_score)
            end
        end
    end
    pulse_alg.optimal_path = pulse_alg.current_optimal_path
    pulse_alg.optimal_objective = pulse_alg.current_optimal_objective
end

function order_nodes(pulse_alg::Pulse, link_dict::Dict{Int, Link}, exploration_order::Function)
    ordered_nodes = sort(collect(keys(link_dict)), by=x->exploration_order(pulse_alg, x)) 
    return ordered_nodes
end

function init_info(deterministic_weights::Vector{String}, 
                   random_weights::Dict{String, Vector{String}})
    deterministic_info = Dict{String, Float64}()
    random_info = Dict{String, Dict{String, Float64}}()
    for weight in deterministic_weights
        deterministic_info[weight] = 0.0
    end
    for random_variable in keys(random_weights)
        random_info[random_variable] = Dict{String, Float64}()
        for weight in random_weights[random_variable]
            random_info[random_variable][weight] = 0.0
        end
    end
    return deterministic_info, random_info
end

function get_path_information(deterministic_weights, random_weights, info_update, path)
    deterministic_info, random_info = init_info(deterministic_weights, random_weights)
    for i in 1:(length(path) - 1)
        deterministic_info, random_info = info_update(path[i], path[i + 1], deterministic_info, random_info)
    end
    return deterministic_info, random_info
end