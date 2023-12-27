# Stochastic pulse class
mutable struct SPulseGraph
    G::Graph
    Primal_Bound::Float64
    Fpath::Vector{Int}
    Resource::Float64
    minimum_time::Int
    T_max::Int
    alpha::Float64
    source::Int
    target::Int
    anim::Vector{Any}
    pos::Vector{Tuple}
    time_expm::Int
    bound::Int
    infeas::Int
    dom::Int
    time_limit::Float64
    pulse_time::Int
    feasibel::Bool
end

function SPulseGraph(G, T_max, alpha, source, target, pos=[])
    return SPulseGraph(G, Inf64, [], 0, 0, T_max, alpha, source, target, [], pos, 0, 0, 0, 0, Inf, 0, false)
end

# Preprocess the graph labeling every node with the shortest cost to the end node (target)
function preprocess(sp::SPulseGraph)
    p = dijkstra_shortest_paths(sp.G, sp.target).dists
    attrs = Dict(i => Dict("labels"=>[], "s_cost"=>p[i]) for i in keys(p))
    attrs.update(Dict(i => Dict("labels"=>[], "s_cost"=>Inf) for i in vertices(sp.G) if i ∉ keys(p)))
    set_props!(sp.G, attrs)

    p = dijkstra_shortest_paths(sp.G, sp.target, edge_type=:min).dists
    attrs = Dict(i => Dict("s_time"=>p[i]) for i in keys(p))
    attrs.update(Dict(i => Dict("s_time"=>Inf) for i in vertices(sp.G) if i ∉ keys(p)))
    set_props!(sp.G, attrs)

    sp.minimum_time = attrs[sp.source]["s_time"]

    for i in vertices(sp.G)
        sp.G[i]["labels"] = []
    end
end

# Check Dominance
function C_Dominance(sp::SPulseGraph, vk, c, prob)
    bool = true
    for (i, (cc, probb)) in enumerate(sp.G[vk]["labels"])
        if c > cc && prob < probb
            bool = false
            sp.dom += 1
        end
    end
    return bool
end

# Check Feasibility
function C_Feasibility(sp::SPulseGraph, vk, prob)
    bool = true
    if prob < sp.alpha
        bool = false
        sp.infeas += 1
    end
    return bool
end

# Check Bounds
function C_Bounds(sp::SPulseGraph, vk, c)
    bool = true
    if c + sp.G[vk]["s_cost"] > sp.Primal_Bound
        bool = false
        sp.bound += 1
    end
    return bool
end

# Update the labels of a given node vk
function update_labels(sp::SPulseGraph, vk, c, prob)
    push!(sp.G[vk]["labels"], (c, prob))
end

function sort(sp::SPulseGraph, sons)
    if sp.feasibel
        return sort(sons, by=x -> sp.G[x]["s_cost"])
    else
        return sort(sons, by=x -> sp.G[x]["s_time"])
    end
end

function update_primal_bound(sp::SPulseGraph, c, prob, P, tRV)
    if prob ≥ sp.alpha && c ≤ sp.Primal_Bound
        sp.Primal_Bound = c
        sp.Fpath = P
        sp.Resource = prob
        sp.tRV = tRV
        sp.feasibel = true
        # println("Primal_Bound: ", c, prob)
    end
end

function pulse(sp::SPulseGraph, vk, c, tRV, tmin, ExptT, P)
    cota = (ExptT - tmin) / (sp.T_max - tmin)  # Makkovitz inequality
    if C_Bounds(sp, vk, c) && vk ∉ P
        if cota < 1 - sp.alpha
            prob = 1 - cota
            # println("Se usa la cota de ", prob)
        else
            t = time()
            prob = cdf(tRV, sp.T_max - tmin - sp.G[vk]["s_time"])
            sp.time_expm += time() - t
            # println("NO se usa la cota de ", 1 - cota)
        end

        update_labels(sp, vk, c, prob)
        if vk == sp.target
            update_primal_bound(sp, c, prob, P + [vk], tRV)
        elseif C_Feasibility(sp, vk, prob) && C_Dominance(sp, vk, c, prob)
            PP = copy(P)
            push!(PP, vk)
            for i in sort(sp, successors(sp.G, vk))
                cc = c + sp.G[vk, i]["Cost"]
                ntRV = tRV + sp.G[vk, i]["tRV"]
                ntmin = tmin + sp.G[vk, i]["tmin"]
                nExptT = ExptT + sp.G[vk, i]["Time"]
                pulse(sp, i, cc, ntRV, ntmin, nExptT, PP)
            end
        end
    end
end

function run_pulse(sp::SPulseGraph, t_limit=1000)
    sp.time_limit = t_limit
    sp.Fpath = []
    sp.Primal_Bound = Inf
    sp.Resource = 0
    preprocess(sp)

    if sp.G[sp.source]["s_cost"] != Inf
        sp.pulse_time = time()
        pulse(sp, vk=sp.source, c=0, tRV=sp.tRV, tmin=0, ExptT=0, P=sp.Fpath)
        sp.pulse_time = time() - sp.pulse_time
    else
        println("The instance is infeasible")
    end

    return sp.Fpath, sp.Primal_Bound, sp.Resource
end

function prob_path(sp::SPulseGraph, path, T_max)
    RV = sp.G[path[1], path[2]]["tRV"]
    arcs = zip(path[2:end-1], path[3:end])

    tmin = 0

    for (i, j) in arcs
        RV += sp.G[i, j]["tRV"]
        tmin += sp.G[i, j]["tmin"]
    end

    return cdf(RV, T_max - tmin)
end


