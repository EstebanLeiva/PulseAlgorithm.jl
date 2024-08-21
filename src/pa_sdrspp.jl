"""
    PaSdrspp

A struct to store the information of the PaSdrspp algorithm.

# Parameters
- `G::Graph`: Graph.
- `α::Float64`: Input reliability.
- `covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}`: Covariance information.
- `variance_costs::Vector{Float64}`: Minimum variance to reach the target node.
- `mean_costs::Vector{Float64}`: Minimum mean to reach the target node.
- `optimal_path::Vector{Int}`: Optimal path.
- `B::Float64`: Best quantile found so far (corresponds to the α-quantile of the optimal path).
- `source_node::Int`: Source node.
- `target_node::Int`: Target node.


"""
mutable struct PaSdrspp
    G::Graph
    α::Float64
    covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}
    variance_costs::Vector{Float64} 
    mean_costs::Vector{Float64} 
    optimal_path::Vector{Int}
    B::Float64 
    source_node::Int 
    target_node::Int
    link_dominance::Dict{Tuple{Int, Int}, PriorityQueue{Tuple{Float64, Float64}, Float64}}
    instance_info::Dict{String, Any}
    pulse_queue::PriorityQueue{Vector{Int}, Float64}
    max_pulse_depth::Int
end

"""
    initialize_PaSdrspp(G::Graph, α::Float64, covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, source_node::String, target_node::String, max_pulse_depth::Int)

Initialize the PaSdrspp struct with a pulse queue with length of exploration max_pulse_depth.
"""
function initialize_PaSdrspp(G::Graph, α::Float64, covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, source_node::String, target_node::String, max_pulse_depth::Int)
    instance_info = Dict(
        "pruned_by_bounds" => 0,
        "total_length_pruned_by_bounds" => 0,
        "pruned_by_dominance" => 0,
        "total_length_pruned_by_dominance" => 0,
        "number_nondominanted_paths" => 0
        )
    source_node = G.name_to_index[source_node]
    target_node = G.name_to_index[target_node]
    return PaSdrspp(G, α, covariance_dict, Vector{Float64}(), Vector{Float64}(), Vector{Int}(), Inf64, source_node, target_node, 
                    PriorityQueue{Tuple{Int, Int}, PriorityQueue{Tuple{Float64, Float64}, Float64}}(), instance_info,
                    PriorityQueue{Vector{Int}, Float64}(), max_pulse_depth)
end

"""
    preprocess!(sdp::PaSdrspp)

Preprocess the graph to calculate the minimum mean and variance to reach the target node.

Note that this information could be saved into a file if the graph is static and could be reused in future runs.
"""
function preprocess!(sdp::PaSdrspp)
    sdp.variance_costs = dijkstra(sdp.G, sdp.target_node, "variance")
    if sdp.variance_costs[sdp.source_node] == Inf
        error("The source node is not reachable from the target node")
    end
    sdp.mean_costs = dijkstra(sdp.G, sdp.target_node, "mean")
end

"""
    check_bounds(sdp::PaSdrspp, current_node::Int, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})

Check if the partial path can improve the current solution. 
    
Return true if it can improve the current solution; and false otherwise.
"""
function check_bounds(sdp::PaSdrspp, current_node::Int, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    mean = mean_path + sdp.mean_costs[current_node]
    variance =  variance_path + covariance_term_path + sdp.variance_costs[current_node]
    dist = Normal(mean, √variance)
    prob = cdf(dist, sdp.B)
    if mean <= sdp.B && prob < sdp.α
        sdp.instance_info["pruned_by_bounds"] += 1
        sdp.instance_info["total_length_pruned_by_bounds"] += length(path)
        return false
    end
    if mean > sdp.B && sdp.α >= 0.5
        sdp.instance_info["pruned_by_bounds"] += 1
        sdp.instance_info["total_length_pruned_by_bounds"] += length(path)
        return false
    end
    quant = quantile(dist, sdp.α)
    if current_node == sdp.target_node && quant <= sdp.B
        sdp.B = quant
        new_path = copy(path)
        push!(new_path, current_node)
        sdp.optimal_path = new_path
        sdp.instance_info["number_nondominanted_paths"] += 1
    end
    return true
end

function check_dominance(sdp::PaSdrspp, current_node::Int, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    if length(path) == 0
        return true
    end
    mean = mean_path
    variance =  variance_path + covariance_term_path
    if (path[end], current_node) ∉ keys(sdp.link_dominance)
        sdp.link_dominance[(path[end], current_node)] = PriorityQueue{Tuple{Float64, Float64}, Float64}()
        enqueue!(sdp.link_dominance[(path[end], current_node)], (mean, variance), mean)
        return true
    end
    queue = sdp.link_dominance[(path[end], current_node)]
    iter_queue = copy(queue)
    while !isempty(iter_queue)
        (mean_iter, variance_iter) = dequeue!(iter_queue)
        if mean_iter < mean && variance_iter < variance
            println(mean_iter, variance_iter)
            println(mean, variance)
            sdp.instance_info["pruned_by_dominance"] += 1
            sdp.instance_info["total_length_pruned_by_dominance"] += length(path)
            return false
        end
    end
    if length(queue) < 10
        enqueue!(queue, (mean, variance), mean)
    end
    return true
end

"""
    pulse(sdp::PaSdrspp, current_node::Int, cost::Float64, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})

Propagate the pulses through the graph while checking the pruning criteria.
"""
function pulse(sdp::PaSdrspp, current_node::Int, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int}, pulse_depth::Int)
    if check_bounds(sdp, current_node, mean_path, variance_path, covariance_term_path, path)
        #if check_dominance(sdp, current_node, mean_path, variance_path, covariance_term_path, path)
            push!(path, current_node)
            link_dict = sdp.G.nodes[current_node].links 
            if path[end] ≠ sdp.target_node
                if pulse_depth < sdp.max_pulse_depth
                    ordered_reachable_nodes = sort(collect(keys(link_dict)), by=x->sdp.mean_costs[x])
                    for reachable_node in ordered_reachable_nodes
                        if reachable_node ∉ path
                            inside_path = copy(path)
                            mean_path_copy = mean_path + link_dict[reachable_node].mean
                            variance_path_copy = variance_path + link_dict[reachable_node].variance
                            covariance_term_path_copy = covariance_term_path + get_covariance_term(sdp.covariance_dict, reachable_node, inside_path)
                            pulse(sdp, reachable_node, mean_path_copy, variance_path_copy, covariance_term_path_copy, inside_path, pulse_depth + 1)
                        end
                    end
                else
                    best_mean = mean_path + sdp.mean_costs[path[end]]
                    best_variance =  variance_path + covariance_term_path + sdp.variance_costs[path[end]]
                    dist = Normal(best_mean, √best_variance)
                    quant = quantile(dist, sdp.α)
                    enqueue!(sdp.pulse_queue, path, quant)
                end
            end
        #end
    end
end

"""
    run_pulse(sdp::PaSdrspp, optimal_path = Vector{Int}(), B = Inf)

Run the PA-SD-RSPP algorithm. Using an initial path and quantile is optional 
and it can be used to improve the performance of the algorithm.

# Arguments
- `sdp::PaSdrspp`: PaSdrspp struct.
- `optimal_path::Vector{Int}`: Initial path.
- `B::Float64`: Best quantile found so far (corresponds to the α-quantile of the initial path).

"""
function run_pulse(sdp::PaSdrspp, optimal_path = Vector{Int}(), B = Inf)
    path = Vector{Int}()
    sdp.optimal_path = optimal_path
    sdp.B = B
    pulse(sdp, sdp.source_node, 0.0, 0.0, 0.0, path, 0)
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