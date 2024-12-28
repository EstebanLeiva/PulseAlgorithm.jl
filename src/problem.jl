struct Problem
    graph::Graph
    source_node::Int
    target_node::Int
    minimization::Bool
    constants::Dict{String, Float64}
end

function Problem(graph::Graph, 
                 source_node::String, 
                 target_node::String, 
                 minimization::Bool, 
                 constants::Dict{String, Float64})
    source_node = graph.name_to_index[source_node]
    target_node = graph.name_to_index[target_node]
    return Problem(graph, source_node, target_node, minimization, constants)
end

function Problem(json_dir::String)
    #TODO: load a problem from json 
    return 1
end