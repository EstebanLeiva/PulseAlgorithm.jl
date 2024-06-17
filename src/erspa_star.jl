mutable struct ERSPAstar
    G::Graph
    α::Float64
    covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64} 
    node_coordinates::Vector{Tuple{Float64, Float64}}
    node_distances::Vector{Float64}
    optimal_path::Vector{Int}
    source_node::Int 
    target_node::Int
    instance_info::Dict{String, Any}
    non_dominated_paths::Dict{Tuple{Int, Int}, PriorityQueue{Vector{Int}, Float64}}
    SE::PriorityQueue{Vector{Int}, Float64}
end

function create_ERSPAstar(G::Graph, α::Float64, covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, source_node::String, target_node::String)
    instance_info = Dict(
        "pruned_by_bounds" => 0,
        "pruned_by_feasibility" => 0, 
        "total_length_pruned_by_bounds" => 0,
        "total_length_pruned_by_feasibility" => 0,
        "number_nondominanted_paths" => 0
        )
    source_node = G.name_to_index[source_node]
    target_node = G.name_to_index[target_node]
    return ERSPAstar(G, α, covariance_dict, Vector{Tuple{Float64, Float64}}(), Vector{Float64}(), Vector{Int}(), source_node, target_node, instance_info, Dict{Tuple{Int, Int}, PriorityQueue{Vector{Int}, Float64}}(), PriorityQueue{Vector{Int}, Float64}())
end

function preprocess_ERSPAstar!()
    return false
end

function euclidean_distance!(erspa::ERSPAstar)
    for (node, _ ) in G.nodes
        erspa.node_distances[node] = sqrt((erspa.node_coordinates[node][1] - erspa.node_coordinates[erspa.target_node][1])^2 + (erspa.node_coordinates[node][2] - erspa.node_coordinates[erspa.target_node][2])^2)
    end
end

function F(erspa::ERSPAstar, path::Vector{Int})
    mean = 0.0
    variance_term = 0.0
    covariance_term = 0.0
    for i in 1:length(path) - 1
        mean += erspa.G.nodes[path[i]].links[path[i+1]].mean
        variance_term += erspa.G.nodes[path[i]].links[path[i+1]].variance
        if length(path)>1 && i < length(path) - 1
            covariance_term += erspa.covariance_dict[(path[i], path[i+1], path[i+1], path[i+2])]
        end
    end
    dist = Normal(mean, √(variance_term + covariance_term))
    f = quantile(dist, erspa.α) + erspa.node_distances[path[end]]
    return f
end

function check_MC_link(erspa::ERSPAstar, link::Tuple{Int, Int} )
    if isempty(keys(erspa.G.nodes[link[2]].links))
        return false
    end
    pos = 0
    neg = 0
    for adjacent in keys(erspa.G.nodes[link[2]].links)
        cov = erspa.covariance_dict[(link[1], link[2], link[2], adjacent)]
        if cov > 0
            pos += 1
        elseif cov < 0
            neg += 1
        end
        if pos > 0 && neg > 0
            return true
        end
    end
    return false     
end

function check_PC_link(erspa::ERSPAstar, link::Tuple{Int, Int})
    if isempty(keys(erspa.G.nodes[link[2]].links))
        return true
    end
    for adjacent in keys(erspa.G.nodes[link[2]].links)
        cov = erspa.covariance_dict[(link[1], link[2], link[2], adjacent)]
        if cov < 0
            return false
        end
    end
    return true
end

function check_NC_link(erspa::ERSPAstar, link::Tuple{Int, Int})
    if isempty(keys(erspa.G.nodes[link[2]].links))
        return false
    end
    for adjacent in keys(erspa.G.nodes[link[2]].links)
        cov = erspa.covariance_dict[(link[1], link[2], link[2], adjacent)]
        if cov > 0
            return false
        end
    end
    return true
end

function check_MB_dominance()
    return false
end

function check_dominance()
    return false
end

function initialization(erspa::ERSPAstar)
    for (reachable_node, link) in erspa.G.nodes[erspa.source_node].links
        path = [erspa.source_node, reachable_node]
        enqueue!(erspa.non_dominated_paths[(erspa.source_node, reachable_node)], path, link.mean)
        enqueue!(erspa.SE, path, F(erspa, path))
    end
end

function path_selection(erspa::ERSPAstar)
    if isempty(erspa.SE)
        return false, []
    end
    path = dequeue!(erspa.SE)
    if path[end] == erspa.target_node
        return false, path
    end
    return true, path
end

function path_extension(erspa::ERSPAstar, path::Vector{Int})
    return false
end