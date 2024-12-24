struct Problem
    graph::Graph, 
    source_node::Int,
    target_node::Int,
    minimization::Bool,
    constants::Dict{String, Float64}
end

function Problem(G::Graph, 
                 source_node::String, 
                 target_node::Stirng, 
                 minimization::Bool, 
                 constants::Dict{String, Float64})
    source_node = G.name_to_index[source_node]
    target_node = G.name_to_index[target_node]
    return Problem(G, source_node, target_node, minimization, constants)
end

function Problem(json_dir::String)
    #TODO: load a problem from json 
    return 1
end