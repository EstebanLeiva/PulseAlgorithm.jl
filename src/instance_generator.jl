function compute_covariance(link1::Link, link2::Link)
    covariance = rand(Uniform(0,1))
    return covariance
end

function pair_links(graph::Graph)
    links1 = get_links(graph)
    links2 = get_links(graph)
    paired_links = Set{Tuple{Int, Int, Int, Int}}()
    count = 0
    for link1 in links1 
        for link2 in links2 
            start_node1, end_node1 = link1
            start_node2, end_node2 = link2
            push!(paired_links, (start_node1, end_node1, start_node2, end_node2))
             #=
            if mod(length(paired_links),100000) == 0 
                create_csv(paired_links, graph, filename)
                paired_links = Set{Tuple{Int, Int, Int, Int}}()
                count += 1
                println("Progress: $count")
            end 
            =#
        end
    end
    #create_csv(paired_links, graph, filename)

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

#=
function run_instance(start_node::String, target_node::String, instance_num::String)
    create_csv(link_pairs, graph, "data/correlation_ChicagoSketch$instance_num.csv")

    dataloader = load_data("data/ChicagoSketch_net.tntp", "ChicagoSketch", "data/correlation_ChicagoSketch$instance_num.csv")

    pulse = create_SPulseGraph(dataloader.graph, 0.9, dataloader.covariance_dict, start_node, target_node, 200.0)
    preprocess(pulse)
    pulse2 =create_SPulseGraph(dataloader.graph, 0.9, dataloader.covariance_dict, start_node, target_node, maximum(pulse.minimum_costs))
    preprocess(pulse2)
    elapsed_time = @elapsed begin
        run_pulse(pulse2)
    end
    
    return elapsed_time, pulse2.instance_info
end
=#

#=
graph = load_graph_from_ta("data/ChicagoSketch_net.tntp","ChicagoSketch")
link_pairs = pair_links(graph)
=#

#=
function load_covariance_dictionary(link_pairs)
    #The CSV file should have the columns start_node_1, end_node_1, start_node_2, end_node_2, covariance
    covariance_dict = Dict{Tuple{Int, Int, Int, Int}, Float64}()

    for (start_node_1, link_index_1, start_node_2, link_index_2) in link_pairs
        link_1 = graph.nodes[start_node_1].links[link_index_1]
        link_2 = graph.nodes[start_node_2].links[link_index_2]
        
        covariance_dict[(start_node_1, link_index_1, start_node_2, link_index_2)] = compute_covariance(link_1, link_2)
    end
    return covariance_dict
end
=#

#=
covariance_dict = load_covariance_dictionary(link_pairs)
dataloader = load_data("data/ChicagoSketch_net.tntp", "ChicagoSketch", covariance_dict)

pulse0 = create_SPulseGraph(dataloader.graph, 0.9, dataloader.covariance_dict, "761", "376", 200.0)
preprocess!(pulse0)
pulse2 =create_SPulseGraph(dataloader.graph, 0.9, dataloader.covariance_dict,  "761", "376", 1*maximum(pulse0.mean_costs))
preprocess!(pulse2)
elapsed_time = @elapsed begin
    run_pulse(pulse2)
end
pulse2.instance_info
=#

function run_instance(start_node::String, target_node::String, instance_num::String)
    dataloader = load_data("data/ChicagoSketch_net.tntp", "ChicagoSketch", "data/correlation_ChicagoSketch$instance_num.csv")

    pulse = create_SPulseGraph(dataloader.graph, 0.9, dataloader.covariance_dict, start_node, target_node, 200.0)
    preprocess(pulse)
    pulse2 =create_SPulseGraph(dataloader.graph, 0.9, dataloader.covariance_dict, start_node, target_node, 3*maximum(pulse.mean_costs))
    preprocess(pulse2)
    elapsed_time = @elapsed begin
        run_pulse(pulse2)
    end
    
    return elapsed_time, pulse2.instance_info
end

#=
results = []
list = [("761","376"), ("217","268"), ("897", "477"), ("274", "84"), 
        ("478", "448"), ("818", "70"), ("918","159"), ("788", "488"), 
        ("663", "902"), ("865", "757"), ("221", "195"), ("517", "444"),
        ("143", "417"), ("542", "750")]
count = 1
for (start, target) in list
    pruned_bound = 0
    pruned_feasibility = 0
    elapsed = 0
    for i in 1:5
        elapsed_time, info = run_instance(start, target, string(count)*"-"*string(i))
        pruned_bound += info["pruned_by_bounds"]
        pruned_feasibility += info["pruned_by_feasibility"]
        elapsed += elapsed_time
    end
    pruned_bound = pruned_bound/5
    pruned_feasibility = pruned_feasibility/5
    elapsed = elapsed/5
    push!(results, (elapsed, pruned_bound, pruned_feasibility))
    println(count)
    count += 1
end

for result in results
    println(result)
end

elapsed_time, info = run_instance("761", "376", "1"*"-"*"1")

for e in results
    println(e)
end
=#

#=
graph = load_graph_from_ta("data/Austin_net.tntp","Austin")

link_pairs = pair_links(graph, "data/correlation_Austin.csv")

dataloader = load_data("data/Austin_net.tntp", "Austin", "data/correlation_Austin.csv")

pulse2 = create_SPulseGraph(dataloader.graph, 0.9, dataloader.covariance_dict, "1", "4218", 200.0)
preprocess(pulse2)
pulse2.mean_costs

elapsed_time = @elapsed begin
    run_pulse(pulse2)
end

pulse2.instance_info
println("Elapsed time: $elapsed_time seconds")
pulse2.B
=#