# Stochastic pulse class
mutable struct SPulseGraph
    G::Graph
    Primal_Bound::Float64
    T_max::Float64
    source_node::String
    target_node::String
end

function SPulseGraph(G::Graph,source_node::String,target_node::String,T_max::Float64)
    return SPulseGraph(G,Inf64,source_node,target_node,T_max)
end

# Preprocess the graph labeling every node with the shortest mean and variance paths to the end node (target)
function preprocess(sp::SPulseGraph)
    #Calculate this shortest_paths values and store them on the dictionary node_labels
end

# Check Dominance
function C_Dominance(sp::SPulseGraph, vk, c, prob)
    
end

# Check Feasibility (true if feasible, false otherwise)
function C_Feasibility(sp::SPulseGraph, node_k, mean_path,variance_path)
    bool = true
    #Replace numbers with shortest mean path and shortes variance path
    dist = Normal(mean_path+1, variance_path+0.1)
    if cdf(dist,sp.T_max) < sp.alpha
        bool = false
    end
    return bool
end

# Check Bounds (true if less than Primal_Bound, false otherwise)
function C_Bounds(sp::SPulseGraph, node_k, cost)
    bool = true
    if c + sp.Optimal_path > sp.Primal_Bound
        bool = false
        sp.bound += 1
    end
    return bool
end

function update_primal_bound(sp::SPulseGraph, c, prob, P, tRV)
    if prob ≥ sp.alpha && c ≤ sp.Primal_Bound
        sp.Primal_Bound = c
        sp.Fpath = P
        sp.Resource = prob
        sp.tRV = tRV
        sp.feasibel = true
        # println("Primal_Bound: ", c, prob)
    end
end

function pulse(sp::SPulseGraph, cost::Float64, mean_path::Float64, variance_path::Float64, path::Vector{Int})
    if C_Feasibility(sp)
end

function run_pulse(sp::SPulseGraph)
    sp.Optimal_path = []
    sp.Primal_Bound = Inf
    preprocess(sp)
    pulse(sp,0,0,0,sp.Optimal_path)
    return sp.Optimal_path, sp.Primal_Bound
end



