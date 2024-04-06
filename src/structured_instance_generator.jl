function modified_dfs2(graph::Graph, start_link::Tuple{Int, Int}, max_depth::Int, depth::Int, visited_pairlinks::Dict{Tuple{Int, Int}, Int}, visited_nodes::Dict{Int, Vector{Int}})
   if depth > max_depth
        return visited_pairlinks
    end

   for adjacent in keys(graph.nodes[start_link[1]].links)
        if !haskey(visited_pairlinks, (start_link[1], adjacent)) || visited_pairlinks[(start_link[1], adjacent)] > depth
            visited_pairlinks[(start_link[1], adjacent)] = depth
            modified_dfs(graph, (start_link[1], adjacent), max_depth, depth + 1, visited_pairlinks, visited_nodes)
        end
    end

    return visited_pairlinks
end

function modified_dfs(graph::Graph, start_link::Tuple{Int, Int}, max_depth::Int, depth::Int, visited_pairlinks::Dict{Tuple{Int, Int}, Int}, previous_node::Int)
    if depth > max_depth
        return visited_pairlinks
    end

    for adjacent in keys(graph.nodes[start_link[2]].links)
        if previous_node != adjacent
            if haskey(visited_pairlinks, (start_link[2], adjacent))
                visited_pairlinks[(start_link[2], adjacent)] = min(visited_pairlinks[(start_link[2], adjacent)], depth)
            else
                visited_pairlinks[(start_link[2], adjacent)] = depth
            end
            modified_dfs(graph, (start_link[2], adjacent), max_depth, depth + 1, visited_pairlinks, adjacent)
        end
    end       

    return visited_pairlinks
end

function get_covariance_dict(graph::Graph, ρ::Float64, max_depth::Int)
    covariance_dict = Dict{Tuple{Int, Int, Int, Int}, Float64}()
    links = get_links_info(graph)
    for link in keys(links)
        visited_pairlinks = modified_dfs(graph, link, max_depth, 1, Dict{Tuple{Int, Int}, Int}(), -1)
        for pairlink in keys(visited_pairlinks)
            if pairlink != link 
                covariance_dict[(link[1], link[2], pairlink[1], pairlink[2])] = (ρ/(visited_pairlinks[pairlink])) * √(links[(pairlink)][3]) * √(links[link][3]) #Corredor structure of ρ's
            end
            if pairlink == link
                covariance_dict[(link[1],link[2],pairlink[1],pairlink[2])] = 1.0
            end
        end
    end

    for link1 in keys(links)
        for link2 in keys(links)
            if !haskey(covariance_dict, (link1[1], link1[2], link2[1], link2[2]))
                covariance_dict[(link1[1], link1[2], link2[1], link2[2])] = 0.0
            end
        end
    end

    return covariance_dict
end

function get_random_pair(graph::Graph)
    nodes = collect(keys(graph.nodes))
    start_node = nodes[rand(1:length(nodes))]
    target_node = nodes[rand(1:length(nodes))]
    return start_node, target_node
end

function get_quantile_path(graph::Graph, path::Vector{Int}, cov_dict::Dict{Tuple{Int, Int, Int, Int}, Float64})
    mean = 0.0
    variance = 0.0
    covariance_term = 0.0

    for i in 1:length(path)-1
        mean += graph.nodes[path[i]].links[path[i+1]].mean
        variance += graph.nodes[path[i]].links[path[i+1]].variance
        for ii in length(path)-1
            if ii != i
                covariance_term += 2*cov_dict[(path[i], path[i+1], path[ii], path[ii+1])]
            end
        end
    end

    covariance_term = 2 * covariance_term
    return mean, variance, covariance_term
end

function get_timeBudget(graph::Graph, start_node::Int, target_node::Int, α::Float64, γ::Float64, cov_dict::Dict{Tuple{Int, Int, Int, Int}, Float64})
    shortest_mean_path = dijkstra_between2Nodes(graph, start_node, target_node, "mean")
    shortest_cost_path = dijkstra_between2Nodes(graph, start_node, target_node, "cost")
    
    mean, variance, covariance_term = get_quantile_path(graph, shortest_mean_path, cov_dict)
    dist = Normal(mean, √(variance+covariance_term))
    T_t_α = quantile(dist, α)

    mean, variance, covariance_term = get_quantile_path(graph, shortest_cost_path, cov_dict)
    dist = Normal(mean, √(variance+covariance_term))
    T_c_α = quantile(dist, α)
    
    T = T_t_α + (T_c_α - T_t_α) * (1 - γ)
    return T
end


function run_structured_instance(start_node::Int, target_node::Int, CV::Float64, ρ::Float64, α::Float64, γ::Float64, max_depth::Int)
    graph = load_graph_from_ta("data/SiouxFalls_net.tntp", "SF", CV)
    covariance_dict = get_covariance_dict(graph, ρ, max_depth)
    T = 2*get_timeBudget(graph, start_node, target_node, α, γ, covariance_dict)
    pulse = create_SPulseGraph(graph, α, covariance_dict, string(start_node),string(target_node),T)
    preprocess!(pulse)
    #save variance cost list as a csv
    CSV.write("variance_costs.csv", DataFrame(variance_costs = pulse.variance_costs), writeheader = false)
    CSV.write("mean_costs.csv", DataFrame(mean_costs = pulse.mean_costs), writeheader = false)
    CSV.write("minimum_costs.csv", DataFrame(minimum_costs = pulse.minimum_costs), writeheader = false)

    elapsed_time = @elapsed begin
        run_pulse(pulse)
    end
    
    return elapsed_time, pulse.instance_info, (start_node, target_node), T, pulse.optimal_path
end

function run_experiment_ρ(start_node::Int, target_node::Int, CV::Float64, ρs::Vector{Float64}, α::Float64, γ::Float64, max_depth::Int)
    results = []
    for ρ in ρs
        elapsed_time, instance_info, pair, T = run_structured_instance(start_node, target_node, CV, ρ, α, γ, max_depth)
        push!(results, (elapsed_time, instance_info, pair, T))
    end
    return results
end

function run_experiment_γ(start_node::Int, target_node::Int, CV::Float64, ρ::Float64, α::Float64, γs::Vector{Float64}, max_depth::Int)
    results = []
    for γ in γs
        elapsed_time, instance_info, pair, T = run_structured_instance(start_node, target_node, CV, ρ, α, γ, max_depth)
        push!(results, (elapsed_time, instance_info, pair, T))
    end
    return results
end

function run_experiment_α(start_node::Int, target_node::Int, CV::Float64, ρ::Float64, αs::Vector{Float64}, γ::Float64, max_depth::Int)
    results = []
    for α in αs
        elapsed_time, instance_info, pair, T = run_structured_instance(start_node, target_node, CV, ρ, α, γ, max_depth)
        push!(results, (elapsed_time, instance_info, pair, T))
    end
    return results
end

function run_experiment_CV(start_node::Int, target_node::Int, CVs::Vector{Float64}, ρ::Float64, α::Float64, γ::Float64, max_depth::Int)
    results = []
    for CV in CVs
        elapsed_time, instance_info, pair, T = run_structured_instance(start_node, target_node, CV, ρ, α, γ, max_depth)
        push!(results, (elapsed_time, instance_info, pair, T))
    end
    return results
end

function run_experiment_max_depth(start_node::Int, target_node::Int, max_depths::Vector{Int}, CV::Float64, ρ::Float64, α::Float64, γ::Float64)
    results = []
    for max_depth in max_depths
        elapsed_time, instance_info, pair, T = run_structured_instance(start_node, target_node, CV, ρ, α, γ, max_depth)
        push!(results, (elapsed_time, instance_info, pair, T))
    end
    return results
end