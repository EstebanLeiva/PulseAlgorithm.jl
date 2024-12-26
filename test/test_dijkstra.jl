using PulseAlgorithm: Graph, Node, create_node!, add_link!, dijkstra, find_node, reverse_graph

@testset "Dijkstra Test (Distances)" begin
    G = Graph(Dict{Int, Node}(), Dict{String, Int}(), Dict{Tuple{String, Tuple{Int, Int, Int, Int}}, Float64}())

    create_node!(G, "1")
    create_node!(G, "2")
    create_node!(G, "3")
    create_node!(G, "4")
    create_node!(G, "5")
    create_node!(G, "s") #6
    create_node!(G, "e") #7

    add_link!(G, "s", "1", Dict("cost" => 2.0), Dict("time" => Dict("mean" => 2.0, "variance" => 3.0)))
    add_link!(G, "1", "e", Dict("cost" => 3.0), Dict("time" => Dict("mean" => 2.0, "variance" => 0.5)))
    add_link!(G, "s", "2", Dict("cost" => 3.0), Dict("time" => Dict("mean" => 2.0, "variance" => 1.0)))
    add_link!(G, "2", "e", Dict("cost" => 5.0), Dict("time" => Dict("mean" => 9.0, "variance" => 1.0)))
    add_link!(G, "s", "3", Dict("cost" => 2.0), Dict("time" => Dict("mean" => 1.0, "variance" => 0.5)))
    add_link!(G, "3", "e", Dict("cost" => 4.0), Dict("time" => Dict("mean" => 1.0, "variance" => 0.5)))
    add_link!(G, "s", "4", Dict("cost" => 1.0), Dict("time" => Dict("mean" => 2.0, "variance" => 3.0)))
    add_link!(G, "4", "5", Dict("cost" => 1.0), Dict("time" => Dict("mean" => 3.0, "variance" => 3.0)))
    add_link!(G, "5", "e", Dict("cost" => 1.0), Dict("time" => Dict("mean" => 2.0, "variance" => 2.0)))

    G_rev = reverse_graph(G)
    node1 = G.name_to_index["s"]
    node2 = G.name_to_index["1"]

    @test G_rev.nodes[node2].links[node1].random["time"]["mean"] == 2.0

    output_dijkstra = dijkstra(G, 7, "cost", false)
    @test output_dijkstra == [3.0, 5.0, 4.0, 2.0, 1.0, 3.0, 0.0]

    output_dijkstra = dijkstra(G, 7, "time" ,"mean", false)
    @test output_dijkstra == [2.0, 9.0, 1.0, 5.0, 2.0, 2.0, 0.0]

    output_dijkstra = dijkstra(G, 7, "time", "variance", false)
    @test output_dijkstra == [0.5, 1.0, 0.5, 5.0, 2.0, 1.0, 0.0]

    _, paths = dijkstra(G, 7, "cost", true)
    @test paths[6] == [6, 4, 5, 7]
end

@testset "Dijkstra Test (Not Connected)" begin
    G = Graph(Dict{Int, Node}(), Dict{String, Int}(), Dict{Tuple{String, Tuple{Int, Int, Int, Int}}, Float64}())

    create_node!(G, "1")
    create_node!(G, "2")
    create_node!(G, "3")
    create_node!(G, "4")
    create_node!(G, "5")
    create_node!(G, "s") #6
    create_node!(G, "e") #7
    create_node!(G, "N") #8

    add_link!(G, "s", "1", Dict("cost" => 2.0), Dict("time" => Dict("mean" => 2.0, "variance" => 3.0)))
    add_link!(G, "1", "e", Dict("cost" => 3.0), Dict("time" => Dict("mean" => 2.0, "variance" => 0.5)))
    add_link!(G, "s", "2", Dict("cost" => 3.0), Dict("time" => Dict("mean" => 2.0, "variance" => 1.0)))
    add_link!(G, "2", "e", Dict("cost" => 5.0), Dict("time" => Dict("mean" => 9.0, "variance" => 1.0)))
    add_link!(G, "s", "3", Dict("cost" => 2.0), Dict("time" => Dict("mean" => 1.0, "variance" => 0.5)))
    add_link!(G, "3", "e", Dict("cost" => 4.0), Dict("time" => Dict("mean" => 1.0, "variance" => 0.5)))
    add_link!(G, "s", "4", Dict("cost" => 1.0), Dict("time" => Dict("mean" => 2.0, "variance" => 3.0)))
    add_link!(G, "4", "5", Dict("cost" => 1.0), Dict("time" => Dict("mean" => 3.0, "variance" => 3.0)))
    add_link!(G, "5", "e", Dict("cost" => 1.0), Dict("time" => Dict("mean" => 2.0, "variance" => 2.0)))


    @test find_node(G, "N") == 8

    output_dijkstra = dijkstra(G, 7, "cost", false)
    @test output_dijkstra == [3.0, 5.0, 4.0, 2.0, 1.0, 3.0, 0.0, Inf]

    output_dijkstra = dijkstra(G, 8, "cost", false)
    @test output_dijkstra == [Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.0]

    _, paths = dijkstra(G, 8, "cost", true)
    @test paths[1] == []
end