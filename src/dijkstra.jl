"""
    dijkstra(graph::Graph, target_node::Int, path::String)

Return a vector with the cost to reach the target node from all other nodes in the graph.

Priority Queue implementation of Dijkstra's algorithm.
"""
function dijkstra(graph::Graph, target_node::Int, path::String)
    n = length(graph.nodes)
    cost = fill(Inf, n)
    cost[target_node] = 0.0
    queue = PriorityQueue()
    enqueue!(queue, target_node, 0.0)
    while !isempty(queue)
        current_node = dequeue!(queue)
        for (neighbor, link) in graph.nodes[current_node].incoming_links
            if path == "variance"
                new_cost = cost[current_node] + link.variance
            elseif path == "mean"
                new_cost = cost[current_node] + link.mean
            elseif path == "cost"
                new_cost = cost[current_node] + link.cost
            else
                error("Unsupported path: $path. Choose 'variance', 'mean', or 'cost'.")
            end
            if new_cost < cost[neighbor]
                cost[neighbor] = new_cost
                if haskey(queue, neighbor)
                    queue[neighbor] = new_cost
                else
                    enqueue!(queue, neighbor, new_cost)
                end
                #queue[neighbor] = queue[neighbor] - new_cost
                #enqueue!(queue, neighbor, new_cost)
            end
        end
    end
    return cost
end

"""
    dijkstra_between_nodes(graph::Graph, start_node::Int, target_node::Int, type::String)

Return the shortest path between two nodes in the graph.

Priority Queue implementation of Dijkstra's algorithm between two nodes.
"""
function dijkstra_between_nodes(graph::Graph, start_node::Int, target_node::Int, type::String)
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
            if type == "mean"
                alt = dist[u] + graph.nodes[u].links[v].mean
            elseif type == "cost"
                alt = dist[u] + graph.nodes[u].links[v].cost
            elseif type == "variance"
                alt = dist[u] + graph.nodes[u].links[v].variance
            else
                error("Unsupported type: $type. Choose 'mean', ''variance', or cost'.")
            end
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