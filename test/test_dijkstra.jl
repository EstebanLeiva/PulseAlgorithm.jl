@testset "Dijkstra Test" begin
    G = Graph(Dict{Int, Node}(), Dict{String, Int}())

    create_node!(G, "1")
    create_node!(G, "2")
    create_node!(G, "3")
    create_node!(G, "4")
    create_node!(G, "5")
    create_node!(G, "s") #6
    create_node!(G, "e") #7

    add_link!(G, "s", "1", 2.0, 2.0, 3.0)
    add_link!(G, "1", "e", 3.0, 2.0, 0.5)
    add_link!(G, "s", "2", 3.0, 2.0, 1.0)
    add_link!(G, "2", "e", 5.0, 9.0, 1.0)
    add_link!(G, "s", "3", 2.0, 1.0, 0.5)
    add_link!(G, "3", "e", 4.0, 1.0, 0.5)
    add_link!(G, "s", "4", 1.0, 2.0, 3.0)
    add_link!(G, "4", "5", 1.0, 3.0, 3.0)
    add_link!(G, "5", "e", 1.0, 2.0, 2.0)

    output_dijkstra = PA.dijkstra(G, 7, "cost")
    @test output_dijkstra == [3.0, 5.0, 4.0, 2.0, 1.0, 3.0, 0.0]

    output_dijkstra = PA.dijkstra(G, 7, "mean")
    @test output_dijkstra == [2.0, 9.0, 1.0, 5.0, 2.0, 2.0, 0.0]

    output_dijkstra = PA.dijkstra(G, 7, "variance")
    @test output_dijkstra == [0.5, 1.0, 0.5, 5.0, 2.0, 1.0, 0.0]

    output_dijkstra_between_nodes = PA.dijkstra_between_nodes(G, 6, 7, "cost")
    @test output_dijkstra_between_nodes == [6, 4, 5, 7]
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

    add_link!(G, "s", "1", 2.0, 2.0, 3.0)
    add_link!(G, "1", "e", 3.0, 2.0, 0.5)
    add_link!(G, "s", "2", 3.0, 2.0, 1.0)
    add_link!(G, "2", "e", 5.0, 9.0, 1.0)
    add_link!(G, "s", "3", 2.0, 1.0, 0.5)
    add_link!(G, "3", "e", 4.0, 1.0, 0.5)
    add_link!(G, "s", "4", 1.0, 2.0, 3.0)
    add_link!(G, "4", "5", 1.0, 3.0, 3.0)
    add_link!(G, "5", "e", 1.0, 2.0, 2.0)

    @test PA.find(G, "N") == 8

    output_dijkstra = PA.dijkstra(G, 7, "cost")
    @test output_dijkstra == [3.0, 5.0, 4.0, 2.0, 1.0, 3.0, 0.0, Inf]

    output_dijkstra = PA.dijkstra(G, 8, "cost")
    @test output_dijkstra == [Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.0]

    output_dijkstra_between_nodes = PA.dijkstra_between_nodes(G, 6, 8, "cost")
    @test output_dijkstra_between_nodes == []
end