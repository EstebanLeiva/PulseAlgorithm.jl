"""
    dijkstra(graph::Graph, target_node::Int, path::String)

Return a vector with the deterministic cost to reach the target node from all other nodes in the graph.

Priority Queue implementation of Dijkstra's algorithm.
"""
function dijkstra(graph::Graph, target_node::Int, det_cost::String, paths::Bool)
    graph = reverse_graph(graph)
    n = length(graph.nodes)
    cost = fill(Inf, n) #TODO: use static array
    cost[target_node] = 0.0
    queue = PriorityQueue() #TODO: use a fast implementation of binary heap
    enqueue!(queue, target_node, 0.0)
    predecessors = fill(-1, n)
    while !isempty(queue)
        current_node = dequeue!(queue)
        for (neighbor, link) in graph.nodes[current_node].links
            new_cost = cost[current_node] + link.deterministic[det_cost]
            if new_cost < cost[neighbor]
                cost[neighbor] = new_cost
                predecessors[neighbor] = current_node
                if haskey(queue, neighbor)
                    queue[neighbor] = new_cost
                else
                    enqueue!(queue, neighbor, new_cost)
                end
            end
        end
    end

    if !paths
        return cost
    end
    
    paths_dict = Dict{Int, Vector{Int}}()
    for node in 1:n
        if cost[node] == Inf
            paths_dict[node] = Vector{Int}()
        else
            path = Vector{Int}()
            current = node
            while current !== -1
                push!(path, current)
                current = predecessors[current]
            end
            paths_dict[node] = path
        end
    end

    return cost, paths_dict
end

"""
    dijkstra(graph::Graph, target_node::Int, rand_var::String, rand_cost::String)

Return a vector with the cost (information associated to rand_var) to reach the target node from all other nodes in the graph.

Priority Queue implementation of Dijkstra's algorithm.
"""
function dijkstra(graph::Graph, target_node::Int, rand_var::String, rand_cost::String, paths::Bool)
    graph = reverse_graph(graph)
    n = length(graph.nodes)
    cost = fill(Inf, n) #TODO: use static array
    cost[target_node] = 0.0
    queue = PriorityQueue() #TODO: use a fast implementation of binary heap
    enqueue!(queue, target_node, 0.0)
    predecessors = fill(-1, n)
    while !isempty(queue)
        current_node = dequeue!(queue)
        for (neighbor, link) in graph.nodes[current_node].links
            new_cost = cost[current_node] + link.random[rand_var][rand_cost]
            if new_cost < cost[neighbor]
                cost[neighbor] = new_cost
                predecessors[neighbor] = current_node
                if haskey(queue, neighbor)
                    queue[neighbor] = new_cost
                else
                    enqueue!(queue, neighbor, new_cost)
                end
            end
        end
    end

    if !paths
        return cost
    end

    paths_dict = Dict{Int, Vector{Int}}()
    for node in 1:n
        if cost[node] == Inf
            paths_dict[node] = Vector{Int}()
        else
            path = Vector{Int}()
            current = node
            while current !== -1
                push!(path, current)
                current = predecessor[current]
            end
            paths_dict[node] = path
        end
    end

    return cost, paths_dict
end