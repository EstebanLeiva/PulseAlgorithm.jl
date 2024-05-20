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
    covariance_dict = DefaultDict{Tuple{Int, Int, Int, Int}, Float64}(0.0) #default value of dic is 0.0
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
    #If the code below is commented, the cov matrix will not have an entry if covariance(i,j) = 0
    #=
    for link1 in keys(links)
        for link2 in keys(links)
            if !haskey(covariance_dict, (link1[1], link1[2], link2[1], link2[2]))
                covariance_dict[(link1[1], link1[2], link2[1], link2[2])] = 0.0
            end
        end
    end
    =#
    return covariance_dict
end

function get_random_pair(graph::Graph)
    nodes = collect(keys(graph.nodes))
    start_node = nodes[rand(1:length(nodes))]
    target_node = nodes[rand(1:length(nodes))]
    return start_node, target_node
end

function get_quantile_path(graph::Graph, path::Vector{Int}, cov_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64})
    mean = 0.0
    variance = 0.0
    covariance_term = 0.0
    for i in 1:length(path)-1
        mean += graph.nodes[path[i]].links[path[i+1]].mean
        variance += graph.nodes[path[i]].links[path[i+1]].variance
        for ii in i + 1:length(path)-1
            covariance_term += 2*cov_dict[(path[i], path[i+1], path[ii], path[ii+1])]
        end
    end
    return mean, variance, covariance_term
end

function get_timeBudget(graph::Graph, start_node::Int, target_node::Int, α::Float64, γ::Float64, cov_dict::DefaultDict{Tuple{Int, Int, Int, Int}, Float64})
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


function run_structured_instance(graph::Graph, start_node::String, target_node::String, ρ::Float64, α::Float64, γ::Float64, max_depth::Int)
    covariance_dict = get_covariance_dict(graph, ρ, max_depth)
    T = 1.2*get_timeBudget(graph, graph.name_to_index[start_node], graph.name_to_index[target_node], α, γ, covariance_dict)
    pulse = create_SPulseGraph(graph, α, covariance_dict, start_node, target_node, T)
    preprocess!(pulse)
    
    #CSV.write("variance_costs_ChicagoRegional_" * target_node * ".csv", DataFrame(variance_costs = pulse.variance_costs), writeheader = false)
    #CSV.write("mean_costs_ChicagoRegional_" * target_node * ".csv", DataFrame(mean_costs = pulse.mean_costs), writeheader = false)
    #CSV.write("minimum_costs_ChicagoRegional_" * target_node * ".csv", DataFrame(minimum_costs = pulse.minimum_costs), writeheader = false)
    
    elapsed_time = @elapsed begin
        run_pulse(pulse)
    end
    
    return elapsed_time, pulse.instance_info, (start_node, target_node), T, pulse.optimal_path
end

function write_shortest_paths(graph::Graph, target_node::String, folder_path::String, network_name::String)
    variance_costs = dijkstra(graph, target_node, "variance")
    file_path = joinpath(folder_path, "variance_costs_" * network_name * "_" * target_node * ".csv")
    CSV.write(file_path, DataFrame(variance_costs = variance_costs), writeheader = false)

    mean_costs = dijkstra(graph, target_node, "mean")
    file_path = joinpath(folder_path, "mean_costs_" * network_name * "_" * target_node * ".csv")
    CSV.write(file_path, DataFrame(mean_costs = mean_costs), writeheader = false)

    minimum_costs = dijkstra(graph, target_node, "cost")
    file_path = joinpath(folder_path, "minimum_costs_" * network_name * "_" * target_node * ".csv")
    CSV.write(file_path, DataFrame(minimum_costs = minimum_costs), writeheader = false)
end

function preprocess_experiments(sp::SPulseGraph, folder_path::String, network_name::String, target_node::String)
    file_path = joinpath(folder_path, "minimum_costs_" * network_name * "_" * target_node * ".csv")
    data = CSV.read(file_path, DataFrame, header = false)
    sp.minimum_costs =  data[:, 1] |> collect

    file_path = joinpath(folder_path, "mean_costs_" * network_name * "_" * target_node * ".csv")
    data = CSV.read(file_path, DataFrame, header = false)
    sp.mean_costs =  data[:, 1] |> collect

    file_path = joinpath(folder_path, "variance_costs_" * network_name * "_" * target_node * ".csv")
    data = CSV.read(file_path, DataFrame, header = false)
    sp.variance_costs =  data[:, 1] |> collect

    nodes = sort(collect(sp.G.nodes), by=x->sp.minimum_costs[x[1]], rev = true)
    possible_start_nodes = filter(x->sp.minimum_costs[x[1]] != Inf, nodes)

    sp.source_node = possible_start_nodes[1][2].name #the node farthest from target

    return sp.source_node
end

function run_experiments_time(graph::Graph, source_node::String, target_node::String, ρ::Float64, α::Float64, γ::Float64, max_depth::Int, folder_path::String, network_name::String)
    covariance_dict = get_covariance_dict(graph, ρ, max_depth)
    pulse = create_SPulseGraph(graph, α, covariance_dict, source_node, target_node, 0.0)
    source_node = preprocess_experiments(pulse, folder_path, network_name, target_node) 
    T = (1)*get_timeBudget(graph, pulse.G.name_to_index[source_node], pulse.G.name_to_index[target_node], α, γ, covariance_dict)
    pulse.T_max = T
    elapsed_time = @elapsed begin
        run_pulse(pulse)
    end
    return elapsed_time, pulse.instance_info, (source_node, target_node), T, pulse.optimal_path
end
