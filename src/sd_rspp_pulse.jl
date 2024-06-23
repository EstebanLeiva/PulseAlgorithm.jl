mutable struct SdrsppPulse
    G::Graph
    α::Float64
    covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}
    minimum_costs::Vector{Float64}
    variance_costs::Vector{Float64} 
    mean_costs::Vector{Float64} 
    optimal_path::Vector{Int}
    B::Float64 
    source_node::String 
    target_node::String
    instance_info::Dict{String, Any}
end

function initialize(G::Graph, α::Float64, covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, source_node::String, target_node::String)
    instance_info = Dict(
        "pruned_by_bounds" => 0,
        "total_length_pruned_by_bounds" => 0,
        "number_nondominanted_paths" => 0
        )
    return SdrsppPulse(G, α, covariance_dict, Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Int}(), Inf64, source_node, target_node, instance_info)
end

function preprocess!(sdp::SdrsppPulse)
    sdp.minimum_costs = dijkstra(sp.G, sp.target_node, "cost")
    if sdp.minimum_costs[sp.G.name_to_index[sp.source_node]] == Inf
        error("The source node is not reachable from the target node")
    end
    sdp.variance_costs = dijkstra(sp.G, sp.target_node, "variance")
    sdp.mean_costs = dijkstra(sp.G, sp.target_node, "mean")
end

function modified_check_bounds(sdp::SdrsppPulse, current_node::Int, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    mean = mean_path + sdp.mean_costs[current_node]
    variance =  variance_path + covariance_term_path + sdp.variance_costs[current_node]
    dist = Normal(mean, √variance)
    quant = quantile(dist, sdp.α)
    if mean < sdp.B && quant > sdp.B
        return false
    end
    return true
end

function pulse(sdp::SdrsppPulse, current_node::Int, cost::Float64, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    if modified_check_bounds(sdp, current_node, mean_path, variance_path, covariance_term_path, path)
        push!(path, current_node)
        link_dict = sp.G.nodes[current_node].links 
        if path[end] ≠ sp.G.name_to_index[sp.target_node]
            ordered_reachable_nodes = sort(collect(keys(link_dict)), by=x->sp.mean_costs[x]) # we explore first the nodes with minimum mean to the end node
            for reachable_node in ordered_reachable_nodes
                if reachable_node ∉ path
                    inside_path = copy(path)
                    cost_copy = cost + link_dict[reachable_node].cost
                    mean_path_copy = mean_path + link_dict[reachable_node].mean
                    variance_path_copy = variance_path + link_dict[reachable_node].variance
                    covariance_term_path_copy = covariance_term_path + calculate_covariance_term(sp, reachable_node, inside_path)
                    pulse(sp, reachable_node, cost_copy, mean_path_copy, variance_path_copy, covariance_term_path_copy, inside_path)
                end
            end
        end
    end
end

function run_pulse(sdp::SdrsppPulse)
    path = Vector{Int}()
    sp.optimal_path = dijkstra_between_nodes(sdp.G, sdp.source_node, sdp.target_node, "mean") 
    mean, variance, covariance_term = get_path_distribution(sdp.G, sdp.optimal_path, sdp.covariance_dict)
    dist = Normal(mean, √(variance + covariance_term))
    sp.B = quantile(dist, sdp.α)
    pulse(sp, sp.G.name_to_index[sp.source_node], 0.0, 0.0, 0.0, 0.0, path)
    return sp.optimal_path, sp.B, sp
end

function calculate_covariance_term(sp::SarPulse, reachable_node::Int, path::Vector{Int})
    n = length(path)
    if n > 1
        last_node = path[end]
        current_node = reachable_node
        covariance_sum = 0.0
        for i in 1:n-1
            covariance_sum += 2 * sp.covariance_dict[(path[i], path[i+1], last_node, current_node)]
        end
        return covariance_sum
    else
        return 0.0
    end
end