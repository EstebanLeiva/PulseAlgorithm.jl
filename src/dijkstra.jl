function dijkstra(graph::Graph, target_node::Int, path::String)
    n = length(graph.nodes)
    # Set the cost of the final node to 0 and the cost of all other nodes to infinity.
    cost = fill(Inf, n)
    cost[target_node] = 0.0

    # Initialize a priority queue.
    # The queue will store the nodes to be visited, 
    #  with the node having the lowest cost having the highest priority.
    queue = PriorityQueue()

    # Add the final node to the priority queue.
    enqueue!(queue, target_node, 0.0)

    # While the queue is not empty, 
    #  remove the node with the highest priority (lowest cost). 
    # For each of its neighbors, 
    #  if the cost of getting to the neighbor via the current node 
    #  is less than the previously known cost, update the cost and 
    #  add the neighbor to the queue.

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
                enqueue!(queue, neighbor, -new_cost) #could be +new_cost
            end
        end
    end

    return cost
end

function dijkstra_between_nodes(graph::Graph, start_node::Int, target_node::Int, type::String)
    # Initialization
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

    
    # Dijkstra's algorithm
    while !isempty(Q)
        u = dequeue!(Q)
        if u == target_node
            break
        end
        for v in keys(graph.nodes[u].links)  #neighbors of u
            if type == "mean"
                alt = dist[u] + graph.nodes[u].links[v].mean
            elseif type == "cost"
                alt = dist[u] + graph.nodes[u].links[v].cost
            else
                error("Unsupported type: $type. Choose 'mean' or 'cost'.")
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
            pushfirst!(S, u)  # Insert u at the beginning of S
            u = prev[u]     # Traverse from target to source
        end
    end
    return S
end