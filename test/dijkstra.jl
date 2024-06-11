using PulseAlgorithm
using Test

@testset "Dijkstra Test" begin
    G = Graph(Dict{Int, Node}(), Dict{String, Int}())

    create_node!(G, "1")
    create_node!(G, "2")
    create_node!(G, "3")
    create_node!(G, "4")
    create_node!(G, "5")
    create_node!(G, "s") #6
    create_node!(G, "e") #7

    add_link!(G, "s", "1", 2.0, 2.0, 2.0)
    add_link!(G, "1", "e", 5.0, 5.0, 5.0)
    add_link!(G, "s", "2", 3.0, 3.0, 3.0)
    add_link!(G, "2", "e", 2.0, 2.0, 2.0)
    add_link!(G, "s", "3", 1.0, 1.0, 1.0)
    add_link!(G, "3", "e", 4.0, 4.0, 4.0)
    add_link!(G, "s", "4", 1.0, 1.0, 1.0)
    add_link!(G, "4", "5", 2.0, 2.0, 2.0)
    add_link!(G, "5", "e", 1.0, 1.0, 1.0)

    output_dijkstra = dijkstra(G, "e", "cost")
    @test output_dijkstra == [5.0, 2.0, 4.0, 3.0, 1.0, 4.0, 0.0]

    output_dijkstra_between2Nodes = dijkstra_between2Nodes(G, 6, 7, "cost")
    @test output_dijkstra_between2Nodes == [6, 4, 5, 7]
end

@testset "Dijkstra Not Connected Test" begin
    G = Graph(Dict{Int, Node}(), Dict{String, Int}())

    create_node!(G, "1")
    create_node!(G, "2")
    create_node!(G, "3")
    create_node!(G, "4")
    create_node!(G, "5")
    create_node!(G, "s") #6
    create_node!(G, "e") #7
    create_node!(G, "N") #8

    add_link!(G, "s", "1", 2.0, 2.0, 2.0)
    add_link!(G, "1", "e", 5.0, 5.0, 5.0)
    add_link!(G, "s", "2", 3.0, 3.0, 3.0)
    add_link!(G, "2", "e", 2.0, 2.0, 2.0)
    add_link!(G, "s", "3", 1.0, 1.0, 1.0)
    add_link!(G, "3", "e", 4.0, 4.0, 4.0)
    add_link!(G, "s", "4", 1.0, 1.0, 1.0)
    add_link!(G, "4", "5", 2.0, 2.0, 2.0)
    add_link!(G, "5", "e", 1.0, 1.0, 1.0)

    @test PulseAlgorithm.find(G, "N") == 8

    output_dijkstra = dijkstra(G, "e", "cost")
    @test output_dijkstra == [5.0, 2.0, 4.0, 3.0, 1.0, 4.0, 0.0, Inf]

    output_dijkstra = dijkstra(G, "N", "cost")
    @test output_dijkstra == [Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.0]

    output_dijkstra_between2Nodes = dijkstra_between2Nodes(G, 6, 8, "cost")
    @test output_dijkstra_between2Nodes == []
end