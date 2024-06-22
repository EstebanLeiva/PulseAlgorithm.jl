using DataStructures
using Distributions

@testset "check_MC_NC_PC Test" begin
    G = Graph(Dict{Int, Node}(), Dict{String, Int}())

    create_node!(G, "1")
    create_node!(G, "2")
    create_node!(G, "3")
    create_node!(G, "4")
    create_node!(G, "5")
    create_node!(G, "6") 
    create_node!(G, "7")
    create_node!(G, "8")

    cost = 1.0
    add_link!(G, "1", "2", cost, 1.0, 1.0)
    add_link!(G, "1", "3", cost, 1.0, 1.0)
    add_link!(G, "1", "4", cost, 1.0, 1.0)
    add_link!(G, "2", "5", cost, 3.0, 1.0)
    add_link!(G, "3", "5", cost, 1.6, 1.0)
    add_link!(G, "4", "5", cost, 1.0, 1.0)
    add_link!(G, "5", "6", cost, 2.0, 1.0)
    add_link!(G, "5", "8", cost, 2.0, 1.0)
    add_link!(G, "6", "8", cost, 2.0, 1.0)
    add_link!(G, "6", "7", cost, 2.0, 1.0)
    add_link!(G, "8", "7", cost, 2.0, 1.0)

    cov_dict = DefaultDict{Tuple{Int, Int, Int, Int}, Float64}(0.0)
    cov_dict[(1, 2, 2, 5)] = 1.0
    cov_dict[(1, 3, 3, 5)] = 1.0
    cov_dict[(1, 4, 4, 5)] = 1.0
    cov_dict[(2, 5, 5, 6)] = 1.0
    cov_dict[(2, 5, 5, 8)] = 1.0
    cov_dict[(3, 5, 5, 6)] = -0.5
    cov_dict[(3, 5, 5, 8)] = -0.5
    cov_dict[(4, 5, 5, 6)] = 0.5
    cov_dict[(4, 5, 5, 8)] = 0.5
    cov_dict[(5, 6, 6, 7)] = 0.5
    cov_dict[(5, 6, 6, 8)] = -1.0
    cov_dict[(5, 8, 8, 7)] = 1.0
    cov_dict[(6, 8, 8, 7)] = -1.0


    erspa = PA.create_ERSPAstar(G, 0.9, cov_dict, Vector{Tuple{Float64, Float64}}(), "1", "7")

    @test PA.check_PC_link(erspa, (1, 2)) == true
    @test PA.check_NC_link(erspa, (1, 2)) == false
    @test PA.check_MC_link(erspa, (1, 2)) == false

    @test PA.check_PC_link(erspa, (1, 3)) == true
    @test PA.check_NC_link(erspa, (1, 3)) == false
    @test PA.check_MC_link(erspa, (1, 3)) == false

    @test PA.check_PC_link(erspa, (1, 4)) == true
    @test PA.check_NC_link(erspa, (1, 4)) == false
    @test PA.check_MC_link(erspa, (1, 4)) == false

    @test PA.check_PC_link(erspa, (2, 5)) == true
    @test PA.check_NC_link(erspa, (2, 5)) == false
    @test PA.check_MC_link(erspa, (2, 5)) == false

    @test PA.check_PC_link(erspa, (3, 5)) == false
    @test PA.check_NC_link(erspa, (3, 5)) == true
    @test PA.check_MC_link(erspa, (3, 5)) == false

    @test PA.check_PC_link(erspa, (4, 5)) == true
    @test PA.check_NC_link(erspa, (4, 5)) == false
    @test PA.check_MC_link(erspa, (4, 5)) == false

    @test PA.check_PC_link(erspa, (5, 6)) == false
    @test PA.check_NC_link(erspa, (5, 6)) == false
    @test PA.check_MC_link(erspa, (5, 6)) == true

    @test PA.check_PC_link(erspa, (5, 8)) == true
    @test PA.check_NC_link(erspa, (5, 8)) == false
    @test PA.check_MC_link(erspa, (5, 8)) == false

    @test PA.check_PC_link(erspa, (6, 7)) == true
    @test PA.check_NC_link(erspa, (6, 7)) == false
    @test PA.check_MC_link(erspa, (6, 7)) == false

    @test PA.check_PC_link(erspa, (8, 7)) == true
    @test PA.check_NC_link(erspa, (8, 7)) == false
    @test PA.check_MC_link(erspa, (8, 7)) == false

end

@testset "Erspa initialization" begin
    G = Graph(Dict{Int, Node}(), Dict{String, Int}())

    create_node!(G, "1")
    create_node!(G, "2")
    create_node!(G, "3")
    create_node!(G, "4")
    create_node!(G, "5")
    create_node!(G, "6") 
    create_node!(G, "7")
    create_node!(G, "8")

    cost = 1.0
    add_link!(G, "1", "2", cost, 1.0, 1.0)
    add_link!(G, "1", "3", cost, 1.0, 1.0)
    add_link!(G, "1", "4", cost, 1.0, 1.0)
    add_link!(G, "2", "5", cost, 3.0, 1.0)
    add_link!(G, "3", "5", cost, 1.6, 1.0)
    add_link!(G, "4", "5", cost, 1.0, 1.0)
    add_link!(G, "5", "6", cost, 2.0, 1.0)
    add_link!(G, "5", "8", cost, 2.0, 1.0)
    add_link!(G, "6", "8", cost, 2.0, 1.0)
    add_link!(G, "6", "7", cost, 2.0, 1.0)
    add_link!(G, "8", "7", cost, 2.0, 1.0)

    cov_dict = DefaultDict{Tuple{Int, Int, Int, Int}, Float64}(0.0)
    cov_dict[(1, 2, 2, 5)] = 1.0
    cov_dict[(1, 3, 3, 5)] = 1.0
    cov_dict[(1, 4, 4, 5)] = 1.0
    cov_dict[(2, 5, 5, 6)] = 1.0
    cov_dict[(2, 5, 5, 8)] = 1.0
    cov_dict[(3, 5, 5, 6)] = -0.5
    cov_dict[(3, 5, 5, 8)] = -0.5
    cov_dict[(4, 5, 5, 6)] = 0.5
    cov_dict[(4, 5, 5, 8)] = 0.5
    cov_dict[(5, 6, 6, 7)] = 0.5
    cov_dict[(5, 6, 6, 8)] = -1.0
    cov_dict[(5, 8, 8, 7)] = 1.0
    cov_dict[(6, 8, 8, 7)] = -1.0

    node_coordinates = Vector{Tuple{Float64, Float64}}()
    push!(node_coordinates, (1.0, 2.0))
    push!(node_coordinates, (2.0, 3.0))
    push!(node_coordinates, (2.0, 2.0))
    push!(node_coordinates, (2.0, 1.0))
    push!(node_coordinates, (3.0, 2.0))
    push!(node_coordinates, (4.0, 2.0))
    push!(node_coordinates, (5.0, 2.0))
    push!(node_coordinates, (4.0, 1.0))

    erspa = PA.create_ERSPAstar(G, 0.9, cov_dict, node_coordinates, "1", "7")

    PA.preprocess_ERSPAstar!(erspa)
    PA.initialization!(erspa)

    @test erspa.non_dominated_paths[(1,2)][[1,2]] == 1.0
    @test erspa.non_dominated_paths[(1,3)][[1,3]] == 1.0
    @test erspa.non_dominated_paths[(1,4)][[1,4]] == 1.0
    @test length(erspa.SE) == 3

    dist = Normal(1.0, 1.0)
    quant = quantile(dist, 0.9)
    elem, priority = peek(erspa.SE) 
    
    @test priority == quant + 3.0
end

@testset "Erspa path selection" begin
    
end

@testset "Erspa path extension" begin
    
end

@testset "Erspa Full Test" begin
    
end