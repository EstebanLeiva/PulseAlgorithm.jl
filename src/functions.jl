### These are the functions that the user would need to specify to run the pulse algorithm for their tailored problem_parameters

# Function that specifies how to update the costs of a path when a new node is added
function update_costs(current_node::Int, reachable_node::Int, path::Vector{Int}, path_det_information::Dict{String, Float64}, path_rand_information::Dict{String, Dict{String, Float64}})
    link_dict = pa.G.nodes[current_node].links 
    path_det_information["cost"] += link_dict[reachable_node].deterministic["cost"]
    path_rand_information["time"]["mean"] += link_dict[reachable_node].random["time"]["mean"]
                        variance_path_copy = variance_path + link_dict[reachable_node].variance
                        covariance_term_path_copy = covariance_term_path + get_covariance_term(pa.covariance_dict, reachable_node, inside_path)
    path_information["time"]["mean"] += link.random["time"]["mean"]
    path_information["cost"] += link.deterministic["cost"]
    return path_information
end
# Funcion that specifies how the objective value is updated when a new node is added

# Functions for each problem specific-pruning rule

