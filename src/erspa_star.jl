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

function create_ERSPAstar(G::Graph, α::Float64, covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, node_coordinates::Vector{Tuple{Float64, Float64}}, source_node::String, target_node::String)
    instance_info = Dict(
        "number_nondominanted_paths" => 0
        )
    source_node = G.name_to_index[source_node]
    target_node = G.name_to_index[target_node]
    return ERSPAstar(G, α, covariance_dict, node_coordinates, Vector{Float64}(undef, length(G.nodes)), Vector{Int}(), source_node, target_node, instance_info, Dict{Tuple{Int, Int}, PriorityQueue{Vector{Int}, Float64}}(), PriorityQueue{Vector{Int}, Float64}())
end

function preprocess_ERSPAstar!(erspa::ERSPAstar)
    euclidean_distance!(erspa)
end

function euclidean_distance!(erspa::ERSPAstar)
    for (node, _ ) in erspa.G.nodes
        erspa.node_distances[node] = sqrt((erspa.node_coordinates[node][1] - erspa.node_coordinates[erspa.target_node][1])^2 + (erspa.node_coordinates[node][2] - erspa.node_coordinates[erspa.target_node][2])^2)
    end
end

function F(erspa::ERSPAstar, path::Vector{Int})
    mean, variance_term, covariance_term = get_path_distribution(erspa.G, path, erspa.covariance_dict)
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

function check_MB_dominance(erspa::ERSPAstar, path::Vector{Int})
    link = (path[end-1], path[end])
    mean_u, variance_term_u, covariance_term_u = get_path_distribution(erspa.G, path, erspa.covariance_dict)
    dist_u = Normal(mean_u, √(variance_term_u + covariance_term_u))
    quant_u = quantile(dist_u, erspa.α)

    P = erspa.non_dominated_paths[link]
    P_copy = copy(P)
    while !isempty(P_copy)
        mean_v, variance_term_v, covariance_term_v = get_path_distribution(erspa.G, dequeue!(P_copy), erspa.covariance_dict)
        if mean_u >= mean_v
            dist_v = Normal(mean_v, √(variance_term_v + covariance_term_v))
            if quant_u > quantile(dist_v, erspa.α)
                return true
            end
        else
            break
        end 
    end  
    return false
end

function check_dominance(erspa::ERSPAstar, path::Vector{Int})
     link = (path[end-1], path[end])
     if check_PC_link(erspa, link)
         β = erspa.α
     elseif check_NC_link(erspa, link) || check_MC_link(erspa, link)
         if erspa.α > 0.5
             β = 0.999
         elseif erspa.α == 0.5
             β = 0.5
         else
             β = 0.001
         end
     end

     P_d = Set{Vector{Int}}()
     v = 1
     mean_u, variance_term_u, covariance_term_u = get_path_distribution(erspa.G, path, erspa.covariance_dict)
     dist_u = Normal(mean_u, √(variance_term_u + covariance_term_u))
     if !haskey(erspa.non_dominated_paths, link)
        erspa.non_dominated_paths[link] = PriorityQueue{Vector{Int}, Float64}()
        enqueue!(erspa.non_dominated_paths[link], path, mean_u)
        return false, P_d
     end

     P = erspa.non_dominated_paths[link]
     P_copy = copy(P)
     while !isempty(P_copy)
        mean_v, variance_term_v, covariance_term_v = get_path_distribution(erspa.G, dequeue!(P_copy), erspa.covariance_dict)
        if mean_u >= mean_v
            dist_v = Normal(mean_v, √(variance_term_v + covariance_term_v))
            if quantile(dist_u, β) > quantile(dist_v, β)
                return true, P_d
            end
        else
            break
        end
    end
    enqueue!(P, path, mean_u)

    higher_priority_paths = get_higher_priority_paths(P, path)
    while !isempty(higher_priority_paths)
        (path_v, v) = dequeue!(higher_priority_paths)
        mean_v, variance_term_v, covariance_term_v = get_path_distribution(erspa.G, path_v, erspa.covariance_dict)
        dist_v = Normal(mean_v, √(variance_term_v + covariance_term_v))
        if quantile(dist_u, β) < quantile(dist_v, β)
            push!(P_d, path_v)
            delete!(P, path_v)
        else
            break
        end
    end

    return false, P_d
end

function initialization!(erspa::ERSPAstar)
    for (reachable_node, link) in erspa.G.nodes[erspa.source_node].links
        path = [erspa.source_node, reachable_node]
        erspa.non_dominated_paths[(erspa.source_node, reachable_node)] = PriorityQueue{Vector{Int}, Float64}()
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
    println("Extended")
    return true, path
end

function path_extension(erspa::ERSPAstar, path::Vector{Int}, input_link::Tuple{Int, Int})
    println("MC: ", check_MC_link(erspa, input_link))
    if check_MC_link(erspa, input_link) 
        MB = check_MB_dominance(erspa, path)
    end
    for (reachable_node, link) in erspa.G.nodes[input_link[2]].links
        if check_MC_link(erspa, input_link)
            if MB && erspa.covariance_dict[input_link[1], input_link[2], input_link[2], reachable_node] >= 0
                continue
            end
        end   
        new_path = copy(path)
        push!(new_path, reachable_node)
        f = F(erspa, new_path)
        dominated, P_d = check_dominance(erspa, new_path)
        println("P_d: ", P_d)
        if !dominated
            enqueue!(erspa.SE, new_path, f)
            for dominated_path in P_d
                delete!(erspa.SE, dominated_path)
            end
        end
    end
end

function run_erspa(erspa::ERSPAstar)
    initialization!(erspa)
    while true
        println(erspa.SE)
        ps = path_selection(erspa)
        if ps[1] == false && !isempty(ps[2])
            erspa.optimal_path = ps[2]
            return erspa
        elseif ps[1] == true
            path_extension(erspa, ps[2], (ps[2][end-1], ps[2][end]))
        else
            return "No path found"
        end
    end
end