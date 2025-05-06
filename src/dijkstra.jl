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
                current = predecessors[current]
            end
            paths_dict[node] = path
        end
    end

    return cost, paths_dict
end

"""
    dijkstra_between_nodes(graph::Graph, start_node::Int, target_node::Int, det_cost::String)

Return the shortest path between two nodes in the graph with respect to det_cost.

Priority Queue implementation of Dijkstra's algorithm between two nodes.
"""
function dijkstra(graph::Graph, start_node::Int, target_node::Int, det_cost::String)
    dist = Vector{Float64}()
    prev = Vector{Int}()
    Q = PriorityQueue()
    for v in sort(collect(keys(graph.nodes)))
        push!(prev, -1)
        if v == start_node
            push!(dist, 0)
        else
            push!(dist, Inf)
        end
        enqueue!(Q, v, dist[v])
    end
    while !isempty(Q)
        u = dequeue!(Q)
        if u == target_node
            break
        end
        for v in keys(graph.nodes[u].links)
            alt = dist[u] + graph.nodes[u].links[v].deterministic[det_cost]
            if alt < dist[v]
                dist[v] = alt
                prev[v] = u
                Q[v] = alt
            end
        end
    end
    S = Vector{Int}()
    u = target_node
    if prev[u] != -1 || u == start_node
        while u != -1
            pushfirst!(S, u)  
            u = prev[u]     
        end
    end
    return S
end

"""
    dijkstra_between_nodes(graph::Graph, start_node::Int, target_node::Int, rand_var::String, rand_cost::String)

Return the shortest path between two nodes in the graph with respect to random cost.

Priority Queue implementation of Dijkstra's algorithm between two nodes.
"""
function dijkstra(graph::Graph, start_node::Int, target_node::Int, rand_var::String, rand_cost::String)
    dist = Vector{Float64}()
    prev = Vector{Int}()
    Q = PriorityQueue()
    for v in sort(collect(keys(graph.nodes)))
        push!(prev, -1)
        if v == start_node
            push!(dist, 0)
        else
            push!(dist, Inf)
        end
        enqueue!(Q, v, dist[v])
    end
    while !isempty(Q)
        u = dequeue!(Q)
        if u == target_node
            break
        end
        for v in keys(graph.nodes[u].links)
            alt = dist[u] + graph.nodes[u].links[v].random[rand_var][rand_cost]
            if alt < dist[v]
                dist[v] = alt
                prev[v] = u
                Q[v] = alt
            end
        end
    end
    S = Vector{Int}()
    u = target_node
    if prev[u] != -1 || u == start_node
        while u != -1
            pushfirst!(S, u)  
            u = prev[u]     
        end
    end
    return S
end