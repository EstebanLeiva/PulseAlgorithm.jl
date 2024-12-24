struct Parameters
    bounds_pruning::Bool,
    feasibility_pruning::Bool,
    dominance_pruning::Bool,
    max_pulse_depth::Int,
    path_completion::Bool,
    exploration_order::String, 
    prep_deterministic_weights::Vector{String},
    prep_random_weights::Dict{String, Vector{String}}
end

function Parameters(bounds_pruning::Bool, 
                    feasibility_pruning::Bool, 
                    dominance_pruning::Bool, 
                    max_pulse_depth::Int, 
                    path_completion::Bool, 
                    exploration_order::String, 
                    prep_deterministic_weights::Vector{String}, 
                    prep_random_weights::Dict{String, Vector{String}})
    if exploration_order ∉ [] #TODO: Add the possible values
        error("The exploration order is not valid. Choose one of the following: ")
    end
    return Parameters(bounds_pruning, feasibility_pruning, dominance_pruning, max_pulse_depth, path_completion, exploration_order, prep_deterministic_weights, prep_random_weights)
end

function Parameters(json_dir)
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
    # pruning strategies
    current_optimal_path::Vector{Int}
    curent_objective::Float64
    const dominance::Dict{Tuple{Int, Int}, PriorityQueue{Tuple{Float64, Float64}, Float64}}
    # acceleration strategies
    const pulse_queue::PriorityQueue{Tuple{Vector{Int}, Dict{String, Float64}, Dict{String, Dict{String, Float64}}}, Float64}
end

function Pulse(problem::Problem, parameters::Parameters)
    det_keys, random_keys = get_links_keys(problem.graph)
    if !(parameters.prep_deterministic_weights ⊆ det_keys)
        error("The deterministic weights for preprocessing are not valid. Choose a subset of the following: $det_keys")
    elseif dict_not_subset(parameters.prep_random_weights, problem.graph)
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
                 Dict{Tuple{Int, Int}, 
                 PriorityQueue{Tuple{Float64, Float64}, Float64}}(), 
                 PriorityQueue{Tuple{Vector{Int}, Dict{String, Float64}, Dict{String, Dict{String, Float64}}}, Float64}())
end

function preprocess!(pulse_alg::Pulse)
    for cost in keys(pulse_alg.parameters.prep_deterministic_weights)
        pulse_alg.prep_deterministic_costs[cost] = dijkstra(pulse_alg.problem.graph, pulse_alg.problem.target_node, cost) #TODO: connect dijkstra
        if pulse_alg.prep_deterministic_costs[cost][pulse_alg.problem.source_node] == Inf
            error("The source node is not reachable from the target node")
        end
        #TODO: Add path completion
    end
    for random_variable in keys(pulse_alg.parameters.prep_random_weights)
        for cost in pulse_alg.parameters.prep_random_weights[random_variable]
            pulse_alg.prep_random_costs[random_variable] = Dict{String, Vector{Float64}}()
            pulse_alg.prep_random_costs[random_variable][cost] = dijkstra(pulse_alg.problem.graph, pulse_alg.problem.target_node, cost) #TODO: connect dijkstra
            if pulse_alg.prep_random_costs[random_variable][cost][pulse_alg.problem.source_node] == Inf
                error("The source node is not reachable from the target node")
            end
        end
        #TODO: Add path completion
    end 
end

function propagate_pulse(pulse_alg::Pulse, 
               current_node::Int,
               deterministic_info::Dict{String, Float64},
               random_info::Dict{String, Dict{String, Float64}},
               current_path::Vector{Int}, 
               current_depth::Int,
               pruning_functions::Vector{Function}, 
               info_update::F1, 
               pulse_score::F2) where {F1<:Function, F2<:Function}
    pass = true
    for pruning_function in pruning_functions
        if pruning_function(pulse_alg, current_node, deterministic_info, random_info, current_path)
            pass = false
            break
        end
    end
    if pass
        push!(current_path, current_node)
        link_dict = pulse_alg.problem.graph.nodes[current_node].links
        if path[end] ≠ pulse_alg.problem.target_node
            if current_depth < pulse_alg.parameters.max_pulse_depth
                ordered_reachable_nodes = order_nodes(link_dict, pulse_alg.parameters.exploration_order)
                for reachable_node in ordered_reachable_nodes
                    if reachable_node ∉ current_path
                        new_path = copy(current_path)
                        new_deterministic_info, new_random_info = info_update(current_node, 
                                                                            reachable_node, 
                                                                            deterministic_info, 
                                                                            random_info) #TODO: check if this is a copy of the info
                        pulse(pulse_alg, reachable_node, new_deterministic_info, new_random_info, new_path, current_depth + 1)
                    end
                end
            else
                score = pulse_score(pulse_alg, current_path, deterministic_info, random_info)
                enqueue!(pulse_alg.pulse_queue, (current_path, deterministic_info, random_info), score)
            end
        end
    end
end

function run_pulse(pulse_alg::Pulse,
                   info_update::F1,
                   pulse_score::F2,
                   init_optimal_path::Vector{Int} = Vector{Int}(), 
                   init_objective::Float64 = Inf) where {F1<:Function, F2<:Function}
    path = Vector{Int}()
    pulse_alg.current_optimal_path = init_optimal_path
    pulse_alg.current_objective = init_objective
    deterministic_info, random_info = init_info(pulse_alg.problem.graph)

    pulse(pulse_alg, pulse_alg.source_node, path, path_information)
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
                pulse(pulse_alg, reachable_node, new_deterministic_info, new_random_info, inside_path, 0)
            end
        end
    end




    while !isempty(sdp.pulse_queue)
        path_to_explore = dequeue!(sdp.pulse_queue)
        mean_path_explore, variance_path_explore, covariance_term_path_explore = get_path_distribution(sdp.G, path_to_explore, sdp.covariance_dict)
        link_dict = sdp.G.nodes[path_to_explore[end]].links 
        
        ordered_reachable_nodes = sort(collect(keys(link_dict)), by=x->sdp.mean_costs[x])
        for reachable_node in ordered_reachable_nodes
            if reachable_node ∉ path_to_explore
                inside_path = copy(path_to_explore)
                mean_path_copy = mean_path_explore + link_dict[reachable_node].mean
                variance_path_copy = variance_path_explore + link_dict[reachable_node].variance
                covariance_term_path_copy = covariance_term_path_explore + get_covariance_term(sdp.covariance_dict, reachable_node, inside_path)
                pulse(sdp, reachable_node, mean_path_copy, variance_path_copy, covariance_term_path_copy, inside_path, 0)
            end
        end
    end
    return sdp.optimal_path, sdp.B, sdp
end

function order_nodes(link_dict, exploration_order)
    #TODO: finish
end

function init_info(graph)
    det_keys, rand_keys = get_links_keys(graph)
    deterministic_info = Dict{String, Float64}()
    random_info = Dict{String, Dict{String, Float64}}()
    for key in det_keys
        deterministic_info[key] = 0.0
    end
    for random_variable in rand_keys
        random_info[random_variable] = Dict{String, Float64}()
        for sub_key in rand_keys[random_variable]
            random_info[random_variable][sub_key] = 0.0
        end
    end
    return deterministic_info, random_info
end

function get_path_information(graph, info_update, path)
    deterministic_info, random_info = init_info(graph)
    for i in 1:(length(path) - 1)
        deterministic_info, random_info = info_update(path[i], path[i + 1], deterministic_info, random_info)
    end
    return deterministic_info, random_info
end