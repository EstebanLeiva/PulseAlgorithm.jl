using PulseAlgorithm


α = 0.9
ρ = 1.0
γ = 0.4
max_depth = 2
CV = 0.5
start_node, target_node = (200,540)

graph = PulseAlgorithm.load_graph_from_ta("data/ChicagoRegional_net.tntp", "data/ChicagoRegional_flow.tntp",  "SF", CV)
covariance_dict = PulseAlgorithm.get_covariance_dict(graph, ρ, max_depth)
elapsed_time, instance_info, (start_node, target_node), T, optimal_path = PulseAlgorithm.run_structured_instance(graph, start_node, target_node, ρ, α, γ, max_depth)
