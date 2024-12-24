"""
    Link(deterministic::Dict{String, Float64}, random::Dict{String, Float64})

Simple link structure with deterministic and random link information.
"""
struct Link
    deterministic::Dict{String, Float64}
    random::Dict{String, Dict{String, Float64}}
end

"""
    Node(name::String, links::Dict{Int, Link})

Simple node structure.
"""
struct Node
    name::String
    links::Dict{Int, Link}
end

"""
    Graph(nodes::Dict{Int, Node}, name_to_index::Dict{String, Int})

Simple graph structure.
"""
struct Graph
    nodes::Dict{Int, Node}
    name_to_index::Dict{String, Int}
    covariance::Dict{Tuple{String, Tuple{Int, Int, Int, Int}}, Float64}
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
function add_link!(graph::Graph, src_name::String, dst_name::String, deterministic::Dict{String, Float64}, random::Dict{String, Float64})
    u = find_or_add!(graph, src_name)
    v = find_or_add!(graph, dst_name)
    graph.nodes[u].links[v] = Link(deterministic, random)
end

"""
    get_links_info(graph::Graph)

Return a dictionary with the link information.
"""
function get_links_info(graph::Graph)
    links = Dict{Tuple{Int, Int}, Tuple{Dict{String, Float64}, Dict{String, Float64}}}()
    for (u, node) in graph.nodes
        for (v,link) in node.links
            links[(u, v)] = (link.deterministic, link.random)
        end
    end
    return links
end

"""
    reverse_graph(graph::Graph)
"""
function reverse_graph(graph::Graph)
    new_graph = Graph(Dict{Int, Node}(), Dict{String, Int}())
    for (u, node) in graph.nodes
        for (v, link) in node.links
            add_link!(new_graph, graph.nodes[v].name, graph.nodes[u].name, link.deterministic, link.random)
        end
    end
    return new_graph
end

function get_link_keys(graph::Graph)
    for (u, node) in graph.nodes
        if !isempty(node.links)
            link = node.links[keys(node.links)[1]]
            det_keys = collect(keys(link.deterministic))
            random_variables = collect(keys(link.random))
            rand_keys = Dict{String, Vector{String}}()
            for random_variable in random_variables
                rand_keys[random_variable] = collect(keys(link.random[random_variable]))
            end
            return det_keys, rand_keys
        end
    end
    return Vector{String}(), Dict{String, Vector{String}}()
end