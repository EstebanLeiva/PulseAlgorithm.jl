using PulseAlgorithm

α = 0.9
ρ = 1.0 #this is fixed for our experiments
γ = 0.4
max_depth = 2
CV = 0.8 #this is fixed for our experiments
start_node, target_node = ("478","448")

graph = PulseAlgorithm.load_graph_from_ta("data/ChicagoSketch_net.tntp", "data/ChicagoSketch_flow.tntp",  "SF", CV)
covariance_dict = PulseAlgorithm.get_covariance_dict(graph, ρ, max_depth)
elapsed_time, instance_info, (start_node, target_node), T, optimal_path = PulseAlgorithm.run_structured_instance(graph, start_node, target_node, ρ, α, γ, max_depth)

###### CONBTROLLED EXPERIMENTS ########
α = 0.6
ρ = 1.0 
γ = 0.0
max_depth = 2

ρ = 1.0 #this is fixed for our experiments
CV = 0.8 #this is fixed for our experiments


## Chicago Sketch
folder_path = "data/shortest_paths/CS"
network_name = "CS"
graph = PulseAlgorithm.load_graph_from_ta("data/ChicagoSketch_net.tntp", "data/ChicagoSketch_flow.tntp",  "SF", CV)
covariance_dict = PulseAlgorithm.get_covariance_dict(graph, ρ, max_depth)
top_nodes = sort(collect(graph.nodes), by = x -> length(x[2].links), rev = true)[1:10]

for (_, node) in top_nodes
    PulseAlgorithm.write_shortest_paths(graph, node.name, "data/shortest_paths/CS", "CS")
end
node = top_nodes[1][2].name
elapsed_time, instance_info, (start_node, target_node), T, optimal_path = PulseAlgorithm.run_experiments_time(graph, node, node, ρ, α, γ, max_depth, folder_path, network_name)

cost = 0.0
for i in 1:length(optimal_path)-1
    cost = cost + graph.nodes[optimal_path[i]].links[optimal_path[i+1]].cost 
end



#Chicago Regional
folder_path = "data/shortest_paths/CR"
network_name = "CR"
graph = PulseAlgorithm.load_graph_from_ta("data/ChicagoRegional_net.tntp", "data/ChicagoRegional_flow.tntp",  "SF", CV)
covariance_dict = PulseAlgorithm.get_covariance_dict(graph, ρ, max_depth)
top_nodes = sort(collect(graph.nodes), by = x -> length(x[2].links), rev = true)[1:10]
for (_, node) in top_nodes
    PulseAlgorithm.write_shortest_paths(graph, node.name, folder_path, network_name)
end

node = top_nodes[4][2].name
elapsed_time, instance_info, (start_node, target_node), T, optimal_path = PulseAlgorithm.run_experiments_time(graph, node, node, ρ, α, γ, max_depth, folder_path, network_name)
length(optimal_path)

##Sioux Falls
α = 0.9
ρ = 1.0 #this is fixed for our experiments
γ = 0.4
max_depth = 2
CV = 0.8 #this is fixed for our experiments
start_node, target_node = ("1","8")
graph = PulseAlgorithm.load_graph_from_ta("data/SiouxFalls_net.tntp", "data/SiouxFalls_flow.tntp",  "SF", CV)
shortest_cost_path = PulseAlgorithm.dijkstra_between2Nodes(graph, graph.name_to_index[start_node], graph.name_to_index[target_node], "cost")
graph.nodes[graph.name_to_index[start_node]].links[2].cost
graph.nodes[2].links
graph.nodes[6].links
cost = 32
graph.nodes[1].links[3].cost + graph.nodes[3].links[4].cost + graph.nodes[4].links[5].cost + graph.nodes[5].links[6].cost + graph.nodes[6].links[8].cost

for i in 1:5
    println(i)
end
