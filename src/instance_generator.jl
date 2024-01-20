function compute_covariance(link1::Link, link2::Link)
    covariance = 0.3
    return covariance
end

function pair_links(graph::Graph)
    links1 = get_links(graph)
    links2 = get_links(graph)
    paired_links = Set{Tuple{Int, Int, Int, Int}}()

    for link1 in links1 
        for link2 in links2 
            start_node1, end_node1 = link1
            start_node2, end_node2 = link2
            push!(paired_links, (start_node1, end_node1, start_node2, end_node2))
        end
    end

    return paired_links
end

function create_csv(link_pairs::Set{Tuple{Int, Int, Int, Int}}, graph::Graph, filename::String)
    header = ["start_node_1", "end_node_1", "start_node_2", "end_node_2", "covariance"]
    rows = []

    for (start_node_1, link_index_1, start_node_2, link_index_2) in link_pairs
        link_1 = graph.nodes[start_node_1].links[link_index_1]
        link_2 = graph.nodes[start_node_2].links[link_index_2]

        covariance = compute_covariance(link_1, link_2)
        row = (start_node_1, link_index_1, start_node_2, link_index_2, covariance)
        push!(rows, row)
    end
    df = DataFrame(rows, header)
    CSV.write(filename, df)
end

#graph = load_graph_from_ta("data/SiouxFalls_net.tntp","SiouxFalls")

#link_pairs = pair_links(graph)

#create_csv(link_pairs, graph, "output.csv")

#dataloader = load_data("data/SiouxFalls_net.tntp", "SiouxFalls", "data/output.csv")

#pulse2 = create_SPulseGraph(dataloader.graph, 0.2, dataloader.covariance_dict, "1", "13", 20.0, 5)

#run_pulse(pulse2)