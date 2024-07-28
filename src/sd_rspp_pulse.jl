mutable struct SdrsppPulse
    G::Graph
    α::Float64
    covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}
    minimum_costs::Vector{Float64}
    variance_costs::Vector{Float64} 
    mean_costs::Vector{Float64} 
    optimal_path::Vector{Int}
    B::Float64 
    source_node::Int 
    target_node::Int
    best_quantile_link::Dict{Tuple{Int, Int}, Float64}
    instance_info::Dict{String, Any}
end

function initialize(G::Graph, α::Float64, covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, source_node::String, target_node::String)
    instance_info = Dict(
        "pruned_by_bounds" => 0,
        "total_length_pruned_by_bounds" => 0,
        "pruned_by_dominance" => 0,
        "total_length_pruned_by_dominance" => 0,
        "number_nondominanted_paths" => 0
        )
    source_node = G.name_to_index[source_node]
    target_node = G.name_to_index[target_node]
    return SdrsppPulse(G, α, covariance_dict, Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Int}(), Inf64, source_node, target_node, Dict{Tuple{Int, Int}, Float64}(), instance_info)
end

function preprocess!(sdp::SdrsppPulse)
    sdp.variance_costs = dijkstra(sdp.G, sdp.target_node, "variance")
    if sdp.variance_costs[sdp.source_node] == Inf
        error("The source node is not reachable from the target node")
    end
    sdp.mean_costs = dijkstra(sdp.G, sdp.target_node, "mean")
end

function check_bounds(sdp::SdrsppPulse, current_node::Int, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    mean = mean_path + sdp.mean_costs[current_node]
    variance =  variance_path + covariance_term_path + sdp.variance_costs[current_node]
    dist = Normal(mean, √variance)
    prob = cdf(dist, sdp.B)
    if mean <= sdp.B && prob < sdp.α - 1e-6
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

function pulse(sdp::SdrsppPulse, current_node::Int, cost::Float64, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    if check_bounds(sdp, current_node, mean_path, variance_path, covariance_term_path, path)
        push!(path, current_node)
        link_dict = sdp.G.nodes[current_node].links 
        if path[end] ≠ sdp.target_node
            ordered_reachable_nodes = sort(collect(keys(link_dict)), by=x->sdp.mean_costs[x])
            for reachable_node in ordered_reachable_nodes
                if reachable_node ∉ path
                    inside_path = copy(path)
                    cost_copy = cost + link_dict[reachable_node].cost
                    mean_path_copy = mean_path + link_dict[reachable_node].mean
                    variance_path_copy = variance_path + link_dict[reachable_node].variance
                    covariance_term_path_copy = covariance_term_path + calculate_covariance_term(sdp.covariance_dict, reachable_node, inside_path)
                    pulse(sdp, reachable_node, cost_copy, mean_path_copy, variance_path_copy, covariance_term_path_copy, inside_path)
                end
            end
        end
    end
end

function run_pulse(sdp::SdrsppPulse, optimal_path = Vector{Int}(), B = Inf)
    path = Vector{Int}()
    sdp.optimal_path = optimal_path
    sdp.B = B
    pulse(sdp, sdp.source_node, 0.0, 0.0, 0.0, 0.0, path)
    return sdp.optimal_path, sdp.B, sdp
end