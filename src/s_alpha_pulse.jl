mutable struct SPulseGraph
    G::Graph
    α::Float64 # reliability threshold
    covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64} # Σ
    minimum_costs::Vector{Float64} # m_c
    variance_costs::Vector{Float64} # m_σ^2
    mean_costs::Vector{Float64} # m_μ
    optimal_path::Vector{Int}
    B::Float64 # Primal Bound
    T_max::Float64 # Time budget
    source_node::String 
    target_node::String
    instance_info::Dict{String, Any}
    check_path::Vector{Int}
end

function create_SPulseGraph(G::Graph, α::Float64, covariance_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64}, source_node::String, target_node::String, T_max::Float64)
    instance_info = Dict(
        "pruned_by_bounds" => 0,
        "pruned_by_feasibility" => 0
        )
    return SPulseGraph(G, α, covariance_dict, Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Int}(), Inf64, T_max, source_node, target_node, instance_info, Vector{Int}())
end

# Preprocess the graph labeling every node with the shortest mean and variance paths to the end node (target)
function preprocess!(sp::SPulseGraph)
    sp.minimum_costs = dijkstra(sp.G, sp.target_node, "cost")
    if sp.minimum_costs[sp.G.name_to_index[sp.source_node]] == Inf
        error("The source node is not reachable from the target node")
    end
    sp.variance_costs = dijkstra(sp.G, sp.target_node, "variance")
    sp.mean_costs = dijkstra(sp.G, sp.target_node, "mean")
end

function preprocess2!(sp::SPulseGraph)
    data = CSV.read("minimum_costs_ChicagoRegional.csv", DataFrame)
    sp.minimum_costs =  data[:, 1] |> collect
    data = CSV.read("mean_costs_ChicagoRegional.csv", DataFrame)
    sp.mean_costs =  data[:, 1] |> collect
    data = CSV.read("variance_costs_ChicagoRegional.csv", DataFrame)
    sp.variance_costs =  data[:, 1] |> collect
end

# Check Feasibility (true if feasible, false otherwise)
function C_Feasibility(sp::SPulseGraph, current_node::Int, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    bool = true
    mean = mean_path + sp.mean_costs[current_node]
    variance =  variance_path + covariance_term_path + sp.variance_costs[current_node]
    dist = Normal(mean, √variance)
    prob = cdf(dist, sp.T_max)
    if sp.T_max >= mean && prob < sp.α
        bool = false
        sp.instance_info["pruned_by_feasibility"] = sp.instance_info["pruned_by_feasibility"] + 1
        check = copy(path)
        check = push!(check, current_node)
        if check == sp.check_path[1:length(path)+1]
            println("Pruned by Feasibility 1 FAIL")
            println(current_node)
            println(path)
            println(mean)
            println(variance_path)
            println(covariance_term_path)
            println(prob)
        end
    elseif sp.T_max < mean && sp.α > 0.5
        bool = false
        sp.instance_info["pruned_by_feasibility"] = sp.instance_info["pruned_by_feasibility"] + 1
        check = copy(path)
        check = push!(check, current_node)
        if check == sp.check_path[1:length(path)+1]
            println("Pruned by Feasibility 2 FAIL")
            println(path)
            println(mean)
            println(sp.mean_costs[current_node])
        end
    end
    return bool
end

# Check Bounds (true if less than B, false otherwise)
function C_Bounds(sp::SPulseGraph, current_node::Int, cost::Float64, path::Vector{Int})
    bool = false
    if cost + sp.minimum_costs[current_node] <= sp.B
        if current_node == sp.G.name_to_index[sp.target_node]
            sp.B = cost
            println("New Primal Bound: $cost")
            new_path = copy(path)
            push!(new_path, current_node)
            sp.optimal_path = new_path
        end
        bool = true
    end
    if !bool
        #if path is an initial subpath of the check path, println
        check = copy(path)
        check = push!(check, current_node)
        if check == sp.check_path[1:length(path)+1]
            println("Pruned by Bounds FAIL: $path")
        end
        sp.instance_info["pruned_by_bounds"] = sp.instance_info["pruned_by_bounds"] + 1
    end
    return bool
end


function pulse(sp::SPulseGraph, current_node::Int, cost::Float64, mean_path::Float64, variance_path::Float64, covariance_term_path::Float64, path::Vector{Int})
    if C_Feasibility(sp, current_node, mean_path, variance_path, covariance_term_path, path)
        if C_Bounds(sp, current_node, cost, path)
            push!(path, current_node)
            link_dict = sp.G.nodes[current_node].links 
            if path[end] ≠ sp.G.name_to_index[sp.target_node]
                ordered_reachable_nodes = sort(collect(keys(link_dict)), by=x->sp.minimum_costs[x]) # we explore first the nodes with minimum cost to the end node
                for reachable_node in ordered_reachable_nodes
                    if reachable_node ∉ path
                        
                        inside_path = copy(path)
                        cost_copy = cost + link_dict[reachable_node].cost
                        mean_path_copy = mean_path + link_dict[reachable_node].mean
                        variance_path_copy = variance_path + link_dict[reachable_node].variance
                        covariance_term_path_copy = covariance_term_path + calculate_covariance_term(sp, reachable_node, inside_path)
                        pulse(sp, reachable_node, cost_copy, mean_path_copy, variance_path_copy, covariance_term_path_copy, inside_path)
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
    println("Running Pulse Algorithm")
    println("Source node: $(sp.source_node)")
    println("Target node: $(sp.target_node)")
    println("Time budget: $(sp.T_max)")
    path = Vector{Int}()
    sp.optimal_path = Vector{Int}() #init optimal path as the minimum mean travel time if it is alpha reliable
    sp.B = Inf
    println(sp.B)
    pulse(sp, sp.G.name_to_index[sp.source_node], 0.0, 0.0, 0.0, 0.0, path)
    return sp.optimal_path, sp.B, sp
end