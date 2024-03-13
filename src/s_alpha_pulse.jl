# Stochastic pulse class
mutable struct SPulseGraph
    G::Graph
    α::Float64 # reliability threshold
    covariance_dict::Dict{Tuple{Int, Int, Int, Int}, Float64} # Σ
    minimum_costs::Vector{Float64} # m_c
    variance_costs::Vector{Float64} # m_σ^2
    mean_costs::Vector{Float64} # m_μ
    optimal_path::Vector{Int}
    B::Float64 # Primal Bound
    T_max::Float64 # Time budget
    source_node::String 
    target_node::String
    instance_info::Dict{String, Any}
end

function create_SPulseGraph(G::Graph, α::Float64, covariance_dict::Dict{Tuple{Int, Int, Int, Int}, Float64}, source_node::String, target_node::String, T_max::Float64)
    instance_info = Dict(
        "pruned_by_bounds" => 0,
        "pruned_by_feasibility" => 0
        )
    return SPulseGraph(G, α, covariance_dict, Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Int}(), Inf64, T_max, source_node, target_node, instance_info)
end

# Preprocess the graph labeling every node with the shortest mean and variance paths to the end node (target)
function preprocess!(sp::SPulseGraph)
    sp.variance_costs = dijkstra(sp.G, sp.target_node, "variance")
    sp.mean_costs = dijkstra(sp.G, sp.target_node, "mean")
    sp.minimum_costs = dijkstra(sp.G, sp.target_node, "cost")
end

# Check Feasibility (true if feasible, false otherwise)
function C_Feasibility(sp::SPulseGraph, current_node::Int, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64)
    bool = true
    mean = mean_path + sp.mean_costs[current_node]
    variance =  variance_path + covariance_term_path + sp.variance_costs[current_node]
    sd = √variance
    println("Approx sd path : $sd")
    println("Approx mean path: $mean")
    dist = Normal(mean, √variance)
    prob = cdf(dist, sp.T_max)
    println("Probability: $prob")
    if prob < sp.α
        bool = false
        sp.instance_info["pruned_by_feasibility"] = sp.instance_info["pruned_by_feasibility"] + 1
    end
    return bool
end

# Check Bounds (true if less than B, false otherwise)
function C_Bounds(sp::SPulseGraph, current_node::Int, cost::Float64, path::Vector{Int})
    bool = false

    min_cost = cost + sp.minimum_costs[current_node]
    println("Approx cost: $min_cost")

    if cost + sp.minimum_costs[current_node] <= sp.B 
        if current_node == sp.G.name_to_index[sp.target_node]
            sp.B = cost
            new_path = copy(path)
            push!(new_path, current_node)
            sp.optimal_path = new_path
        end
        bool = true
    end
    if !bool
        sp.instance_info["pruned_by_bounds"] = sp.instance_info["pruned_by_bounds"] + 1
    end
    return bool
end


function pulse(sp::SPulseGraph, current_node::Int, cost::Float64, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    if C_Feasibility(sp, current_node, mean_path, variance_path, covariance_term_path)
        if C_Bounds(sp, current_node, cost, path)
            push!(path, current_node)
            pathend = path[end]
            println("Path END: $pathend")
            link_dict = sp.G.nodes[current_node].links 
            if path[end] ≠ sp.G.name_to_index[sp.target_node]
                
                for reachable_node in keys(link_dict)
                    if reachable_node ∉ path
                        #prints the current node and the path
                        reachable_node_str = string(reachable_node)
                        inside_path_str = join(path, " -> ")
                        println("Reached node $reachable_node_str with path $inside_path_str")
                        
                        inside_path = copy(path)
                        cost = cost + link_dict[reachable_node].cost
                        println("Cost path: $cost")
                        mean_path = mean_path + link_dict[reachable_node].mean
                        println("Mean path: $mean_path")
                        variance_path = variance_path + link_dict[reachable_node].variance
                        println("Variance path: $variance_path")
                        covariance_term_path = covariance_term_path + calculate_covariance_term(sp, reachable_node, inside_path)
                        println("Covariance term path: $covariance_term_path")
                        println(" ")
                        pulse(sp, reachable_node, cost, mean_path, variance_path, covariance_term_path, inside_path)
                    end
                end
            end
        end
    end
end

function calculate_covariance_term(sp::SPulseGraph, reachable_node::Int, path::Vector{Int})
    n = length(path)
    if n > 1
        last_node = path[end]
        current_node = reachable_node
        covariance_sum = 0.0
        for i in 1:n-1
            covariance_sum += 2 * sp.covariance_dict[(path[i], path[i+1], last_node, current_node)]
        end
        return covariance_sum
    else
        return 0.0
    end
end

function run_pulse(sp::SPulseGraph)
    path = Vector{Int}()
    sp.optimal_path = Vector{Int}()
    sp.B = Inf
    pulse(sp, sp.G.name_to_index[sp.source_node], 0.0, 0.0, 0.0, 0.0, path)
    return sp.optimal_path, sp.B, sp
end