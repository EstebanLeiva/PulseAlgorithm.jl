struct Link
    cost::Int
    value::Int
end

mutable struct Node
    name::String
    links::Dict{Int, Link}
end

mutable struct Graph
    nodes::Dict{Int, Node}
    name_to_index::Dict{String, Int}
end

function create_node!(graph::Graph, name::String)
    n = length(graph.nodes) + 1
    node = Node(name, Dict{Int, Link}())
    graph.nodes[n] = node
    graph.name_to_index[name] = n
    return n
end

function find(graph::Graph, name::String)
    return get(graph.name_to_index, name,-1)
end

function find_or_add!(graph::Graph, name::String)
    n = find(graph, name)
    if n < 0
        n = create_node!(graph, name)
    end
    return n
end

function add_link!(graph::Graph, src_name::String, dst_name::String, cost::Float64)
    u = find_or_add!(graph, src_name)
    v = find_or_add!(graph, dst_name)
    graph.nodes[u].links[v] = Link(cost, 0)
    if false # Set to your condition for directed (false) or undirected (true) graph
        graph.nodes[v].links[u] = Link(cost, 0)
    end
end