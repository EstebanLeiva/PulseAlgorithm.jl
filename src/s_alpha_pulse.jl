# Stochastic pulse class
mutable struct SPulseGraph
    G::Graph
    alpha::Float64
    covariance_dict::Dict{Tuple{Int, Int, Int, Int}, Float64}
    variance_costs::Vector{Float64}
    mean_costs::Vector{Float64}
    optimal_path::Vector{Int}
    Primal_Bound::Float64
    T_max::Float64
    source_node::String
    target_node::String
    dominance_set::Dict{Int, Vector{Tuple{Float64, Float64}}}
    dominance_set_size::Int
end

function create_SPulseGraph(G::Graph, alpha::Float64, covariance_dict::Dict{Tuple{Int, Int, Int, Int}, Float64}, source_node::String, target_node::String, T_max::Float64, dominance_set_size::Int)
    dominance_set = Dict{Int, Vector{Tuple{Float64, Float64}}}()
    for node in keys(G.nodes)
        dominance_set[node] = Vector{Tuple{Float64, Float64}}(undef, dominance_set_size)
    end
    return SPulseGraph(G, alpha, covariance_dict, Vector{Float64}(), Vector{Float64}(), Vector{Int}(), Inf64, T_max, source_node, target_node,
                        dominance_set, dominance_set_size)
end

# Preprocess the graph labeling every node with the shortest mean and variance paths to the end node (target)
function preprocess(sp::SPulseGraph)
    sp.variance_costs = dijkstra(sp.G, sp.target_node, "variance")
    sp.mean_costs = dijkstra(sp.G, sp.target_node, "mean")
end

# Check Dominance (true if it is not dominated, false otherwise)
function C_Dominance(sp::SPulseGraph, current_node::Int, cost::Float64, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64)
    dist = Normal(mean_path, variance_path + covariance_term_path)
    prob = cdf(dist, sp.T_max)
    tuple = (cost,prob)
    index = searchsortedfirst(sp.dominance_set[current_node], tuple, by=x->x[1])
    if index < sp.dominance_set_size
        if length(sp.dominance_set[current_node]) < sp.dominance_set_size
            insert!(sp.dominance_set[current_node], index, tuple)
            return true
        else
            for i in index:sp.dominance_set_size
                if prob < sp.dominance_set[current_node][i][2]
                    sp.dominance_set[current_node][i] = tuple
                    sort!(sp.dominance_set[current_node], by = x -> x[1])        
                    return true
                end
            end
            return true
        end
    else
        #we just discard strongly dominated paths
        if prob < sp.dominance_set[current_node][end][2]
            return false
        else
            return true
        end
    end
    
end

# Check Feasibility (true if feasible, false otherwise)
function C_Feasibility(sp::SPulseGraph, current_node::Int, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64)
    bool = true
    dist = Normal(mean_path + sp.mean_costs[current_node], variance_path + covariance_term_path + sp.variance_costs[current_node])
    prob = cdf(dist, sp.T_max)
    if prob < sp.alpha
        bool = false
        println("Pruned by Infeasibility: "*string(prob))
    end
    return bool
end

# Check Bounds (true if less than Primal_Bound, false otherwise)
function C_Bounds(sp::SPulseGraph, current_node::Int, cost::Float64, path::Vector{Int})
    bool = false
    if cost < sp.Primal_Bound
        if current_node == sp.G.name_to_index[sp.target_node]
            sp.Primal_Bound = cost
            new_path = copy(path)
            push!(new_path, current_node)
            sp.optimal_path = new_path
        end
        bool = true
    end
    return bool
end


function pulse(sp::SPulseGraph, current_node::Int, cost::Float64, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    #the probability bound can be calculated here and reused in C_Feasibility and C_Dominance
    if C_Feasibility(sp, current_node, mean_path, variance_path, covariance_term_path)
        if C_Bounds(sp, current_node, cost, path)
            if C_Dominance(sp, current_node, cost, mean_path, variance_path, covariance_term_path)
                push!(path, current_node)
                link_dict = sp.G.nodes[current_node].links 
                if path[end] != sp.G.name_to_index[sp.target_node]
                    println("path: "*string(path))
                    println("Reachable Nodes: "*string(keys(link_dict)))
                    for reachable_node in keys(link_dict)
                        inside_path = copy(path)
                        println("reachable_node: "*string(reachable_node))
                        cost = cost + link_dict[reachable_node].cost
                        println(cost)
                        mean_path = mean_path + link_dict[reachable_node].mean
                        variance_path = variance_path + link_dict[reachable_node].variance
                        covariance_term_path = covariance_term_path + calculate_covariance_term(sp, reachable_node, inside_path)
                        pulse(sp, reachable_node, cost, mean_path, variance_path, covariance_term_path, inside_path)
                    end
                end
            end
        end
    end
end

function calculate_covariance_term(sp::SPulseGraph, reachable_node::Int, path::Vector{Int})
    #Check for error in calculating the covariance for the case of same node
    n = length(path)
    if n > 1
        last_node = path[end]
        current_node = reachable_node
        covariance_sum = 0.0
        for i in 1:n-1
            covariance_sum += sp.covariance_dict[(path[i], path[i+1], last_node, current_node)]
        end
        return 2*covariance_sum
    else
        return 0.0
    end
end

function run_pulse(sp::SPulseGraph)
    path = Vector{Int}()
    sp.optimal_path = Vector{Int}()
    sp.Primal_Bound = Inf
    preprocess(sp)
    pulse(sp, sp.G.name_to_index[sp.source_node], 0.0, 0.0, 0.0, 0.0, path)
    return sp.optimal_path, sp.Primal_Bound
end