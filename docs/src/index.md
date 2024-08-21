# PulseAlgorithm.jl
```@contents
```

## Introduction
PulseAlgorithm.jl provides an implementation for the Pulse Algorithm for Reliable Shortest Path Problems (PA-RSPPs) designed to tackle challenging stochastic shortest path problems in networks with normally distributed link travel times with correlation. This algorithm scales efficiently in large networks and performs competitively with state-of-the-art methods. This package is designed for researchers and practitioners in transportation engineering, logistics, and telecommunications in which these types of networks may arise. 

The package is available under MIT License. 

PulseAlgorithm documentation is organized in two sections:
* Graph Structure: documentation for the graph data structure used.
* Data Loading: documentation for the data loader for real transportation networks found in [Transportation Networks](https://github.com/bstabler/TransportationNetworks)
* Util: documentation for utility functions.
* Pulse Algorithm for Reliable Shortest Path Problems (PA-RSPPs): documentation for the PA-RSPPs algorithm applied to different settings.


## Installation

To install PulseAlgorithm.jl, you can use the Julia package manager. Open the Julia REPL and run the following command:
```
pkg> add https://github.com/EstebanLeiva/PulseAlgorithm.jl
```

## Graph

```@docs
PA.Link
PA.Node
PA.Graph
PA.create_node!
PA.find
PA.find_or_add!
PA.add_link!
PA.get_links_info
```

### Dijkstra

```@docs
PA.dijkstra
PA.dijkstra_between_nodes
```

## Data Loader
```@docs
PA.TAData
PA.load_ta
PA.load_flowCost
PA.calculate_fft_coefficient
PA.load_graph_from_ta
```
## Util

```@docs
PA.get_path_distribution
PA.modified_dfs
PA.get_covariance_dict
PA.get_higher_priority_paths
PA.get_covariance_term
```

## PA-S-``\alpha``RP

```@docs
PA.PaSarp
PA.initialize_PaSarp
PA.preprocess!(::PaSarp)
PA.check_feasibility(::PaSarp, ::Int, ::Float64, ::Float64, ::Float64, ::Vector{Int})
PA.check_bounds(::PaSarp, ::Int, ::Float64, ::Vector{Int})
PA.pulse(::PaSarp, ::Int, ::Float64, ::Float64, ::Float64, ::Float64, ::Vector{Int})
PA.run_pulse(::PaSarp, ::Vector{Int}, B::Float64)
```

## PA-SD_RSPP
```@docs
PA.PaSdrspp
PA.initialize_PaSdrspp
PA.preprocess!(::PaSdrspp)
PA.check_bounds(::PaSdrspp, ::Int, ::Float64, ::Float64, ::Float64, ::Vector{Int})
PA.pulse(::PaSdrspp, ::Int, ::Float64, ::Float64, ::Float64, ::Vector{Int}, ::Int)
PA.run_pulse(::PaSdrspp, ::Vector{Int}, B::Float64)
```