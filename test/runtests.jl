using PulseAlgorithm
using DataFrames
using CSV

α = 0.9
ρ = 1.0
γ = 0.4
max_depth = 2
CV = 0.8
start_node, target_node = (1,11)

graph = PulseAlgorithm.load_graph_from_ta("data/SiouxFalls_net.tntp", "SF", CV)
covariance_dict = PulseAlgorithm.get_covariance_dict(graph, ρ, max_depth)
elapsed_time, instance_info, (start_node, target_node), T, optimal_path = PulseAlgorithm.run_structured_instance(start_node, target_node, CV, ρ, α, γ, max_depth)

covariance_dict[(1,2,2,6)]
dict = PulseAlgorithm.get_links_info(graph)
#change keys to strings
dict = Dict(string(k) => v for (k,v) in dict)
df = DataFrame(Key = collect(keys(dict)), Value = collect(values(dict)))
#create a new column for each of the entries of the values of the dict
for i in 1:length(df[1,2])
    df[!,Symbol("Link $(i)")] = [df[j,2][i] for j in 1:size(df,1)]
end
CSV.write("output.csv", df)

#=
α = 0.9
ρs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
γ = 0.4
max_depth = 5
CV = 0.8
start_node, target_node = (735,382)

results_ρs = PulseAlgorithm.run_experiment_ρ(start_node, target_node, CV, ρs, α, γ, max_depth)
times_ρs = [results_ρs[i][1] for i in eachindex(results_ρs)]
pruned_bounds_ρs = [results_ρs[i][2]["pruned_by_bounds"] for i in eachindex(results_ρs)]
pruned_feasibility_ρs = [results_ρs[i][2]["pruned_by_feasibility"] for i in eachindex(results_ρs)]

plot(ρs, times_ρs, xlabel="ρ", ylabel="Time (s)", title="Time vs ρ: α=0.9, γ=0.4, radius=5, CV=0.8, Budget = 3*T", lw=2)
plot(ρs, pruned_bounds_ρs, xlabel="ρ", ylabel="Pruned by bounds", title="Bounds vs ρ: α=0.9, γ=0.4, radius=5, CV=0.8, Budget = 3*T", lw=2)
plot(ρs, pruned_feasibility_ρs, xlabel="ρ", ylabel="Pruned by feasibility", title="Feasibility vs ρ: α=0.9, γ=0.4, radius=5, CV=0.8, Budget = 3*T", lw=2)

α = 0.9
γs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
ρ = 0.5
max_depth = 20.0
CV = 0.8
start_node, target_node = (735,382)

results_γs = PulseAlgorithm.run_experiment_γ(start_node, target_node, CV, ρ, α, γs, max_depth)
times_γs = [results_γs[i][1] for i in eachindex(results_γs)]
pruned_bounds_γs = [results_γs[i][2]["pruned_by_bounds"] for i in eachindex(results_γs)]
pruned_feasibility_γs = [results_γs[i][2]["pruned_by_feasibility"] for i in eachindex(results_γs)]

plot(γs, times_γs, xlabel="γ", ylabel="Time (s)", title="Time vs γ: α=0.9, ρ=0.5, radius=20, CV=0.8", lw=2)
plot(γs, pruned_bounds_γs, xlabel="γ", ylabel="Pruned by bounds", title="Pruned by bounds vs γ", lw=2)
plot(γs, pruned_feasibility_γs, xlabel="γ", ylabel="Pruned by feasibility", title="Pruned by feasibility vs γ", lw=2)

αs = [0.6, 0.7, 0.8, 0.9, 0.95]
γ = 0.4
ρ = 0.5
max_depth = 5
CV = 0.8
start_node, target_node = (735,382)

results_αs = PulseAlgorithm.run_experiment_α(start_node, target_node, CV, ρ, αs, γ, max_depth)
times_αs = [results_αs[i][1] for i in eachindex(results_αs)]
pruned_bounds_αs = [results_αs[i][2]["pruned_by_bounds"] for i in eachindex(results_αs)]
pruned_feasibility_αs = [results_αs[i][2]["pruned_by_feasibility"] for i in eachindex(results_αs)]

plot(αs, times_αs, xlabel="α", ylabel="Time (s)", title="Time vs α: γ=0.4, ρ=0.5, radius=5, CV=0.8, Budget = 3T", lw=2)
plot(αs, pruned_bounds_αs, xlabel="α", ylabel="Pruned by bounds", title="Pruned by bounds vs α: γ=0.4, ρ=0.5, radius=5, CV=0.8, Budget = 3T", lw=2)
plot(αs, pruned_feasibility_αs, xlabel="α", ylabel="Pruned by feasibility", title="Pruned by feasibility vs α: γ=0.4, ρ=0.5, radius=5, CV=0.8, Budget = 3T", lw=2)

CVs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.0, 5.0]
ρ = 0.5
α = 0.9
γ = 0.4
max_depth = 5
start_node, target_node = (735,382)

results_CVs = PulseAlgorithm.run_experiment_CV(start_node, target_node, CVs, ρ, α, γ, max_depth)
times_CVs = [results_CVs[i][1] for i in eachindex(results_CVs)]
pruned_bounds_CVs = [results_CVs[i][2]["pruned_by_bounds"] for i in eachindex(results_CVs)]
pruned_feasibility_CVs = [results_CVs[i][2]["pruned_by_feasibility"] for i in eachindex(results_CVs)]

plot(CVs, times_CVs, xlabel="CV", ylabel="Time (s)", title="Time vs CV: α=0.9, γ=0.4, ρ=0.5, radius=5", lw=2)
plot(CVs, pruned_bounds_CVs, xlabel="CV", ylabel="Pruned by bounds", title="Pruned by bounds vs CV", lw=2)
plot(CVs, pruned_feasibility_CVs, xlabel="CV", ylabel="Pruned by feasibility", title="Pruned by feasibility vs CV", lw=2)

max_depths = [1, 2, 3, 4, 5]
CV = 0.8
ρ = 0.5
α = 0.9
γ = 0.4
start_node, target_node = (735,382)

results_max_depths = PulseAlgorithm.run_experiment_max_depth(start_node, target_node, max_depths, CV, ρ, α, γ)
times_max_depths = [results_max_depths[i][1] for i in eachindex(results_max_depths)]
pruned_bounds_max_depths = [results_max_depths[i][2]["pruned_by_bounds"] for i in eachindex(results_max_depths)]
pruned_feasibility_max_depths = [results_max_depths[i][2]["pruned_by_feasibility"] for i in eachindex(results_max_depths)]

plot(max_depths, times_max_depths, xlabel="Max depth", ylabel="Time (s)", title="Time vs Max depth: α=0.9, γ=0.4, ρ=0.5, CV=0.8", lw=2)
plot(max_depths, pruned_bounds_max_depths, xlabel="Max depth", ylabel="Pruned by bounds", title="Pruned by bounds vs Max depth", lw=2)
plot(max_depths, pruned_feasibility_max_depths, xlabel="Max depth", ylabel="Pruned by feasibility", title="Pruned by feasibility vs Max depth", lw=2)
=#