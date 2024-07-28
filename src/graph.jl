"""
    Link(cost::Float64, mean::Float64, variance::Float64, value::Float64)

Simple link structure.

# Parameters
- `cost::Float64`: Link cost.
- `mean::Float64`: Link travel time mean.
- `variance::Float64`: Link travel time variance.
- `value::Float64`: Link value.
"""
mutable struct Link
    cost::Float64
    mean::Float64
    variance::Float64
    value::Float64
end

"""
    Node(name::String, links::Dict{Int, Link}, incoming_links::Dict{Int, Link})

Simple node structure.

The links dictionary stores the outgoing links where the key is the destination node index and the value is the link itself.
The incoming_links dictionary stores the incoming links where the key is the source node index and the value is the link itself.
"""
mutable struct Node
    name::String
    links::Dict{Int, Link}
    incoming_links::Dict{Int, Link}
end

"""
    Graph(nodes::Dict{Int, Node}, name_to_index::Dict{String, Int})

Simple graph structure. 

The nodes are stored in a dictionary where the key is the node index and the value is the node itself. 
The name_to_index dictionary is used to map the node name to the node index.
"""
mutable struct Graph
    nodes::Dict{Int, Node}
    name_to_index::Dict{String, Int}
end

"""
    create_node!(graph::Graph, name::String)

Create a new node in the graph.
"""
function create_node!(graph::Graph, name::String)
    n = length(graph.nodes) + 1
    node = Node(name, Dict{Int, Link}(), Dict{Int, Link}())
    graph.nodes[n] = node
    graph.name_to_index[name] = n
    return n
end

"""
    find(graph::Graph, name::String)

Find the node index given the node name.
"""
function find(graph::Graph, name::String)
    return get(graph.name_to_index, name, -1)
end

"""
    find_or_add!(graph::Graph, name::String)

Find the node index given the node name. If the node does not exist, create a new node.
"""
function find_or_add!(graph::Graph, name::String)
    n = find(graph, name)
    if n < 0
        n = create_node!(graph, name)
    end
    return n
end

"""
    add_link!(graph::Graph, src_name::String, dst_name::String, cost::Float64, mean::Float64, variance::Float64)

Add a link between two nodes.

"""
function add_link!(graph::Graph, src_name::String, dst_name::String, cost::Float64, mean::Float64, variance::Float64)
    u = find_or_add!(graph, src_name)
    v = find_or_add!(graph, dst_name)
    graph.nodes[u].links[v] = Link(cost, mean, variance, 0)
    graph.nodes[v].incoming_links[u] = Link(cost, mean, variance, 0)
end

"""
    get_links_info(graph::Graph)

Return a dictionary with the link information.
"""
function get_links_info(graph::Graph)
    links = Dict{Tuple{Int, Int}, Tuple{Float64, Float64, Float64}}()
    for (u, node) in graph.nodes
        for (v,link) in node.links
            links[(u, v)] = (link.cost, link.mean, link.variance)
        end
    end
    return links
end