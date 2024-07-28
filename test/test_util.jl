@testset "Covariance Dict Test" begin
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

    output_covariance_dict = PA.get_covariance_dict(G, 1.0, 2)
    @test output_covariance_dict[(6, 2, 2, 7)] == 1.0*√(3.0)*√(2.0)
    @test output_covariance_dict[(6, 1, 1, 7)] == 1.0*√(2.0)*√(5.0)
    @test output_covariance_dict[(6, 3, 3, 7)] == 1.0*√(1.0)*√(4.0)
    @test output_covariance_dict[(6, 4, 4, 5)] == 1.0*√(1.0)*√(2.0)
    @test output_covariance_dict[(6, 4, 5, 7)] == (1.0/2)*√(1.0)*√(1.0)
end