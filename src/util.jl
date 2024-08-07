"""
    get_path_distribution(graph::Graph, path::Vector{Int}, cov_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64})

Return the mean, variance and covariance term of a given path's travel time distribution.

See also [`get_covariance_term`](@ref).
"""
function get_path_distribution(graph::Graph, path::Vector{Int}, cov_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64})
    mean = 0.0
    variance = 0.0
    covariance_term = 0.0
    for i in 1:length(path)-1
        mean += graph.nodes[path[i]].links[path[i+1]].mean
        variance += graph.nodes[path[i]].links[path[i+1]].variance
        for ii in i + 1:length(path)-1
            covariance_term += 2*cov_dict[(path[i], path[i+1], path[ii], path[ii+1])]
        end
    end
    return mean, variance, covariance_term
end

"""
    modified_dfs(graph::Graph, start_link::Tuple{Int, Int}, max_depth::Int, depth::Int, visited_pairlinks::Dict{Tuple{Int, Int}, Int}, previous_node::Int)

Return the links at a distance less or equal to max_depth from the starting link.

See also [`get_covariance_dict`](@ref).
"""
function modified_dfs(graph::Graph, start_link::Tuple{Int, Int}, max_depth::Int, depth::Int, visited_pairlinks::Dict{Tuple{Int, Int}, Int}, previous_node::Int)
    if depth > max_depth
        return visited_pairlinks
    end

    for adjacent in keys(graph.nodes[start_link[2]].links)
        if previous_node != adjacent
            if haskey(visited_pairlinks, (start_link[2], adjacent))
                visited_pairlinks[(start_link[2], adjacent)] = min(visited_pairlinks[(start_link[2], adjacent)], depth)
            else
                visited_pairlinks[(start_link[2], adjacent)] = depth
            end
            modified_dfs(graph, (start_link[2], adjacent), max_depth, depth + 1, visited_pairlinks, adjacent)
        end
    end       

    return visited_pairlinks
end

"""
    get_covariance_dict(graph::Graph, ρ::Float64, max_depth::Int)

Return the covariance dictionary of a graph following a spatial correlation structure of radius=max_depth.

See also [`modified_dfs`](@ref).
"""
function get_covariance_dict(graph::Graph, ρ::Float64, max_depth::Int)
    covariance_dict = DefaultDict{Tuple{Int, Int, Int, Int}, Float64}(0.0) #default value of dic is 0.0
    links = get_links_info(graph)
    for link in keys(links)
        visited_pairlinks = modified_dfs(graph, link, max_depth, 1, Dict{Tuple{Int, Int}, Int}(), -1)
        for pairlink in keys(visited_pairlinks)
            if pairlink != link 
                covariance_dict[(link[1], link[2], pairlink[1], pairlink[2])] = (ρ/(visited_pairlinks[pairlink])) * √(links[(pairlink)][3]) * √(links[link][3]) #Corredor structure of ρ's
            end
            if pairlink == link
                covariance_dict[(link[1],link[2],pairlink[1],pairlink[2])] = 1.0
            end
        end
    end
    return covariance_dict
end

"""
    get_higher_priority_paths(pq::PriorityQueue, element::Vector{Int})

Return the elements with higher priority than the given element in a priority queue.
"""
function get_higher_priority_paths(pq::PriorityQueue, element::Vector{Int})
    temp_pq = copy(pq) 
    while !isempty(temp_pq)
        (k, v) = dequeue!(temp_pq)
        if k == element
            return temp_pq
        end
    end
    return temp_pq
end

"""
    get_covariance_term(covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, reachable_node::Int, path::Vector{Int})

Return the covariance term of the given path's travel time distribution.

See also [`get_path_distribution`](@ref), [`pulse`](@ref).
"""
function get_covariance_term(covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, reachable_node::Int, path::Vector{Int})
    n = length(path)
    if n > 1
        last_node = path[end]
        current_node = reachable_node
        covariance_sum = 0.0
        for i in 1:n-1
            covariance_sum += 2 * covariance_dict[(path[i], path[i+1], last_node, current_node)]
        end
        return covariance_sum
    else
        return 0.0
    end
end