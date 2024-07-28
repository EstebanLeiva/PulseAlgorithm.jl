"""
    PaSarp(G::Graph, α::Float64, covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, minimum_costs::Vector{Float64}, variance_costs::Vector{Float64}, mean_costs::Vector{Float64}, optimal_path::Vector{Int}, B::Float64, T_max::Float64, source_node::Int, target_node::Int, instance_info::Dict{String, Any})

Define the struct that will contain the information of the PA-S-αRP algorithm.

# Parameters
- `G::Graph`: Graph.
- `α::Float64`: Reliability constraint.
- `covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}`: Covariance information.
- `minimum_costs::Vector{Float64}`: Minimum cost to reach the target node.
- `variance_costs::Vector{Float64}`: Minimum variance to reach the target node.
- `mean_costs::Vector{Float64}`: Minimum mean to reach the target node.
- `optimal_path::Vector{Int}`: Optimal path.
- `B::Float64`: Best cost found so far (corresponds to the cost of the optimal path).
- `T_max::Float64`: Travel time budget.
- `source_node::Int`: Source node.
- `target_node::Int`: Target node.
- `instance_info::Dict{String, Any}`: Instance information.

"""
mutable struct PaSarp
    G::Graph
    α::Float64
    covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}
    minimum_costs::Vector{Float64}
    variance_costs::Vector{Float64} 
    mean_costs::Vector{Float64} 
    optimal_path::Vector{Int}
    B::Float64 
    T_max::Float64 
    source_node::Int 
    target_node::Int
    instance_info::Dict{String, Any}
end

"""
    initialize_PaSarp(G::Graph, α::Float64, covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, source_node::String, target_node::String, T_max::Float64)

Initialize the PaSarp struct.

# Arguments

- `G::Graph`: Graph.
- `α::Float64`: Reliability constraint.
- `covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}`: Covariance information.
- `source_node::String`: Source node.
- `target_node::String`: Target node.
- `T_max::Float64`: Travel time budget.

"""
function initialize_PaSarp(G::Graph, α::Float64, covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, source_node::String, target_node::String, T_max::Float64)
    instance_info = Dict(
        "pruned_by_bounds" => 0,
        "pruned_by_feasibility" => 0, 
        "total_length_pruned_by_bounds" => 0,
        "total_length_pruned_by_feasibility" => 0,
        "number_nondominanted_paths" => 0
        )
    source_node = G.name_to_index[source_node]
    target_node = G.name_to_index[target_node]
    return PaSarp(G, α, covariance_dict, Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Int}(), Inf64, T_max, source_node, target_node, instance_info)
end

"""
    preprocess!(pa::PaSarp)

Calculate the minimum costs, variance costs, and mean costs to reach the target node; 
while ensuring that there is a path between the source and target nodes.

Note that this information could be saved into a file if the graph is static and could be reused in future runs.
"""
function preprocess!(pa::PaSarp)
    pa.minimum_costs = dijkstra(pa.G, pa.target_node, "cost")
    if pa.minimum_costs[pa.source_node] == Inf
        error("The source node is not reachable from the target node")
    end
    pa.variance_costs = dijkstra(pa.G, pa.target_node, "variance")
    pa.mean_costs = dijkstra(pa.G, pa.target_node, "mean")
end

"""
    check_feasibility(pa::PaSarp, current_node::Int, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})

Check if the path is feasible with respect to the time budget and the reliability constraint. 

Return true if the path is feasible; and false otherwise
"""
function check_feasibility(pa::PaSarp, current_node::Int, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    bool = true
    mean = mean_path + pa.mean_costs[current_node]
    variance =  variance_path + covariance_term_path + pa.variance_costs[current_node]
    dist = Normal(mean, √variance)
    prob = cdf(dist, pa.T_max)
    if pa.T_max >= mean && prob < pa.α
        bool = false
        pa.instance_info["pruned_by_feasibility"] += 1
        pa.instance_info["total_length_pruned_by_feasibility"] += length(path)
        
    elseif pa.T_max < mean && pa.α > 0.5
        bool = false
        pa.instance_info["pruned_by_feasibility"] += 1
        pa.instance_info["total_length_pruned_by_feasibility"] += length(path)
        end
    return bool
end

"""
    check_bounds(pa::PaSarp, current_node::Int, cost::Float64, path::Vector{Int})

Check if the partial path can improve the current solution. 
    
Return true if it can improve the current solution; and false otherwise.
"""
function check_bounds(pa::PaSarp, current_node::Int, cost::Float64, path::Vector{Int})
    bool = false
    if cost + pa.minimum_costs[current_node] <= pa.B
        if current_node == pa.target_node
            pa.B = cost
            new_path = copy(path)
            push!(new_path, current_node)
            pa.optimal_path = new_path
            pa.instance_info["number_nondominanted_paths"] += 1
        end
        bool = true
    end
    if !bool
        pa.instance_info["pruned_by_bounds"] += 1
        pa.instance_info["total_length_pruned_by_bounds"] += length(path)
    end
    return bool
end

"""
    pulse(pa::PaSarp, current_node::Int, cost::Float64, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})

Propagate the pulses through the graph while checking the pruning criteria.
"""
function pulse(pa::PaSarp, current_node::Int, cost::Float64, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    if check_feasibility(pa, current_node, mean_path, variance_path, covariance_term_path, path)
        if check_bounds(pa, current_node, cost, path)
            push!(path, current_node)
            link_dict = pa.G.nodes[current_node].links 
            if path[end] ≠ pa.target_node
                if pa.B != Inf
                    ordered_reachable_nodes = sort(collect(keys(link_dict)), by=x->pa.minimum_costs[x]) # we explore first the nodes with minimum cost to the end node
                else
                    ordered_reachable_nodes = sort(collect(keys(link_dict)), by=x->pa.mean_costs[x]) # we explore first the nodes with minimum mean to the end node
                end
                for reachable_node in ordered_reachable_nodes
                    if reachable_node ∉ path
                        inside_path = copy(path)
                        cost_copy = cost + link_dict[reachable_node].cost
                        mean_path_copy = mean_path + link_dict[reachable_node].mean
                        variance_path_copy = variance_path + link_dict[reachable_node].variance
                        covariance_term_path_copy = covariance_term_path + get_covariance_term(pa.covariance_dict, reachable_node, inside_path)
                        pulse(pa, reachable_node, cost_copy, mean_path_copy, variance_path_copy, covariance_term_path_copy, inside_path)
                    end
                end
            end
        end
    end
end

"""
    run_pulse(pa::PaSarp, optimal_path::Vector{Int}, B::Float64)

Run the PA-S-αRP algorithm. Using an initial optimal path and cost is optional 
and it can be used to improve the performance of the algorithm.

# Arguments
- `pa::PaSarp`: PaSarp struct.
- `optimal_path::Vector{Int}=Vector{Int}()`: Initial path.
- `B::Float64=Inf`: Cost of the initial  path.
"""
function run_pulse(pa::PaSarp, optimal_path::Vector{Int} = Vector{Int}(), B::Float64 = Inf)
    path = Vector{Int}()
    pa.optimal_path = optimal_path #init optimal path as  user-specified if it is alpha reliable
    pa.B = B #init B as the cost of a user-specified path if it is alpha reliable
    pulse(pa, pa.source_node, 0.0, 0.0, 0.0, 0.0, path)
    return pa.optimal_path, pa.B, pa
end