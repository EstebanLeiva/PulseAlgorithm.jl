mutable struct SarPulse
    G::Graph
    α::Float64
    covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}
    minimum_costs::Vector{Float64}
    variance_costs::Vector{Float64} 
    mean_costs::Vector{Float64} 
    optimal_path::Vector{Int}
    B::Float64 
    T_max::Float64 
    source_node::Int 
    target_node::Int
    instance_info::Dict{String, Any}
end

function initialize(G::Graph, α::Float64, covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, source_node::String, target_node::String, T_max::Float64)
    instance_info = Dict(
        "pruned_by_bounds" => 0,
        "pruned_by_feasibility" => 0, 
        "total_length_pruned_by_bounds" => 0,
        "total_length_pruned_by_feasibility" => 0,
        "number_nondominanted_paths" => 0
        )
    source_node = G.name_to_index[source_node]
    target_node = G.name_to_index[target_node]
    return SarPulse(G, α, covariance_dict, Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Int}(), Inf64, T_max, source_node, target_node, instance_info)
end

function preprocess!(sp::SarPulse)
    sp.minimum_costs = dijkstra(sp.G, sp.target_node, "cost")
    if sp.minimum_costs[sp.source_node] == Inf
        error("The source node is not reachable from the target node")
    end
    sp.variance_costs = dijkstra(sp.G, sp.target_node, "variance")
    sp.mean_costs = dijkstra(sp.G, sp.target_node, "mean")
end

function check_feasibility(sp::SarPulse, current_node::Int, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    bool = true
    mean = mean_path + sp.mean_costs[current_node]
    variance =  variance_path + covariance_term_path + sp.variance_costs[current_node]
    dist = Normal(mean, √variance)
    prob = cdf(dist, sp.T_max)
    if sp.T_max >= mean && prob < sp.α
        bool = false
        sp.instance_info["pruned_by_feasibility"] += 1
        sp.instance_info["total_length_pruned_by_feasibility"] += length(path)
        
    elseif sp.T_max < mean && sp.α > 0.5
        bool = false
        sp.instance_info["pruned_by_feasibility"] += 1
        sp.instance_info["total_length_pruned_by_feasibility"] += length(path)
        end
    return bool
end

function check_bounds(sp::SarPulse, current_node::Int, cost::Float64, path::Vector{Int})
    bool = false
    if cost + sp.minimum_costs[current_node] <= sp.B
        if current_node == sp.target_node
            sp.B = cost
            new_path = copy(path)
            push!(new_path, current_node)
            sp.optimal_path = new_path
            sp.instance_info["number_nondominanted_paths"] += 1
        end
        bool = true
    end
    if !bool
        sp.instance_info["pruned_by_bounds"] += 1
        sp.instance_info["total_length_pruned_by_bounds"] += length(path)
    end
    return bool
end

function pulse(sp::SarPulse, current_node::Int, cost::Float64, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    if check_feasibility(sp, current_node, mean_path, variance_path, covariance_term_path, path)
        if check_bounds(sp, current_node, cost, path)
            push!(path, current_node)
            link_dict = sp.G.nodes[current_node].links 
            if path[end] ≠ sp.target_node
                if sp.B != Inf
                    ordered_reachable_nodes = sort(collect(keys(link_dict)), by=x->sp.minimum_costs[x]) # we explore first the nodes with minimum cost to the end node
                else
                    ordered_reachable_nodes = sort(collect(keys(link_dict)), by=x->sp.mean_costs[x]) # we explore first the nodes with minimum mean to the end node
                end
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
end

function run_pulse(sp::SarPulse, optimal_path = Vector{Int}(), B = Inf)
    path = Vector{Int}()
    sp.optimal_path = optimal_path #init optimal path as  user-specified if it is alpha reliable
    sp.B = B #init B as the cost of a user-specified path if it is alpha reliable
    pulse(sp, sp.source_node, 0.0, 0.0, 0.0, 0.0, path)
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