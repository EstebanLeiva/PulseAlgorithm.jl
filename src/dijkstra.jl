function dijkstra(graph::Graph, target_node::String, path::String)
    target_node = graph.name_to_index[target_node]
    n = length(graph.nodes)
    # Set the cost of the final node to 0 and the cost of all other nodes to infinity.
    cost = fill(Inf, n)
    cost[target_node] = 0

    # Initialize a priority queue.
    # The queue will store the nodes to be visited, 
    #  with the node having the lowest cost having the highest priority.
    queue = PriorityQueue()

    # Add the final node to the priority queue.
    enqueue!(queue, target_node, 0)

    # While the queue is not empty, 
    #  remove the node with the highest priority (lowest cost). 
    # For each of its neighbors, 
    #  if the cost of getting to the neighbor via the current node 
    #  is less than the previously known cost, update the cost and 
    #  add the neighbor to the queue.

    while !isempty(queue)
        current_node = dequeue!(queue)

        for (neighbor, link) in graph.nodes[current_node].links
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
                enqueue!(queue, neighbor, -new_cost)
            end
        end
    end

    return cost
end