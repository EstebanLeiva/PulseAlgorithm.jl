# Stochastic pulse class
mutable struct SPulseGraph
    G::Graph
    optimal_path::Vector{Int}
    Primal_Bound::Float64
    T_max::Float64
    source_node::String
    target_node::String
    dominance_set::Dict{Int, Vector{Tuple{Float64, Float64}}}
    dominance_set_size::Int
end

function SPulseGraph(G::Graph,source_node::String,target_node::String,T_max::Float64,dominance_set_size::Int)
    return SPulseGraph(G,Vector{Int}(),Inf64,source_node,target_node,T_max,Dict{Int, Vector{Tuple{Float64, Float64}}}(),dominance_set_size)
end

# Preprocess the graph labeling every node with the shortest mean and variance paths to the end node (target)
function preprocess(sp::SPulseGraph)
    #Calculate this shortest_paths values and store them on the dictionary node_labels
end

# Check Dominance (true if it is not dominated, false otherwise)
function C_Dominance(sp::SPulseGraph, current_node::Int, cost, probability)
    #We need to calculate the full variance + covariances for each path for using dominance?? I think so
    #An idea would be to add the covariance terms as the sum of the greatest covariance of each path
    return true
end

# Check Feasibility (true if feasible, false otherwise)
function C_Feasibility(sp::SPulseGraph, current_node::Int, mean_path,variance_path)
    bool = true
    #change the numbers to the shortest mean and variance paths calculated with the mst
    dist = Normal(mean_path+1, variance_path+0.1)
    if cdf(dist,sp.T_max) < sp.alpha
        bool = false
    end
    return bool
end

# Check Bounds (true if less than Primal_Bound, false otherwise)
function C_Bounds(sp::SPulseGraph, current_node::Int, cost::Float64, path::Vector{Int})
    bool = false
    if cost < sp.Primal_Bound
        if current_node == sp.target_node
            sp.Primal_Bound = cost
            new_path = copy(path)
            push!(new_path, current_node)
            sp.optimal_path = new_path
        end
        bool = true
    end
    return bool
end


function pulse(sp::SPulseGraph, current_node::Int, cost::Float64, mean_path::Float64, variance_path::Float64, path::Vector{Int})
    #the probability bound can be calculated here and reused in C_Feasibility and C_Dominance
    if C_Feasibility(sp, current_node, mean_path, variance_path)
        if C_Bounds(sp, current_node, cost, path)
            if C_Dominance(sp, current_node, cost, probability)
                push!(path, current_node)
                link_dict = sp.G.nodes[current_node].links 
                for reachable_node in keys(link_dict)
                    cost = cost + link_dict[reachable_node].cost
                    mean_path = mean_path + link_dict[reachable_node].mean
                    variance_path = variance_path + link_dict[reachable_node].variance
                    pulse(sp, reachable_node, cost, mean_path, variance_path, path)
                end
            end
        end
    end
end

function run_pulse(sp::SPulseGraph)
    path = []
    sp.optimal_path = []
    sp.Primal_Bound = Inf
    preprocess(sp)
    pulse(sp, sp.source_node, 0, 0, 0, path)
    return sp.optimal_path, sp.Primal_Bound
end



