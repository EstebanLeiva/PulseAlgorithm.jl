using PulseAlgorithm: load_graph_from_ta, find_node
@testset "load_graph_from_ta Test" begin
    CV = 0.1
    toll_factor = 1.0
    distance_factor = 1.0

    graph = PA.load_graph_from_ta("data/SiouxFalls_net.tntp", "data/SiouxFalls_flow.tntp", "SF", CV, toll_factor, distance_factor)

    @test length(graph.nodes) == 24

    @test graph.nodes[1].name == "1"
    @test graph.nodes[2].name == "2"
    @test haskey(graph.nodes[1].links, 2)
    @test haskey(graph.nodes[PA.find_node(graph, "24")].links, PA.find_node(graph,"23"))
    @test graph.nodes[1].links[2].random["time"]["mean"] - 6 * ( 1 + 0.15 * (6.0008162373543197/25900.20064)^4 ) < 1e-3
    @test graph.nodes[1].links[2].random["time"]["variance"] - CV * 6 * ( 1 + 0.15 * (6.0008162373543197/25900.20064)^4 ) < 1e-3
    @test graph.nodes[1].links[2].deterministic["cost"] - (6 * ( 1 + 0.15 * (6.0008162373543197/25900.20064)^4 ) + toll_factor * 0 + distance_factor * 6) < 1e-3
end