graph = load_graph_from_ta("data/SiouxFalls_net.tntp","ChicagoSketch")


function modified_dfs(graph::Graph, start_link::Tuple{Int, Int}, max_depth::Float64, depth::Float64, visited_pairlinks::Dict{Tuple{Int, Int}, Float64})
   if depth >= max_depth
        return visited_pairlinks
    end

   visited_pairlinks[start_link] = 0.0
   for adjacent in keys(graph.nodes[start_link[1]].links)
        if !haskey(visited_pairlinks, (start_link[1], adjacent))
            visited_pairlinks[(start_link[1], adjacent)] = depth
            modified_dfs(graph, (start_link[1], adjacent), max_depth, depth + 1.0, visited_pairlinks)
        end
    end

    return visited_pairlinks
end




function get_covariance_dict(graph::Graph, ρ::Float64)
    covariance_dict = Dict{Tuple{Int, Int, Int, Int}, Float64}()
    links = get_links_info(graph)
    for link in keys(links)
        visited_pairlinks = modified_dfs(graph, link, 3.0, 1.0, Dict{Tuple{Int, Int}, Float64}())
        for pairlink in keys(visited_pairlinks)
            if pairlink != link && visited_pairlinks[pairlink] != 0.0
                covariance_dict[(link[1],link[2],pairlink[1],pairlink[2])] = ρ/(visited_pairlinks[pairlink]) * √(links[(pairlink)][3]) * √(links[link][3])
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
                covariance_dict[(link2[1], link2[2], link1[1], link1[2])] = 0.0
            end
        end
    end

    return covariance_dict
end

#cov_dict = get_covariance_dict(graph, 0.5)
