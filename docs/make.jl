using Documenter, PulseAlgorithm
const PA = PulseAlgorithm

makedocs(
    sitename = "PulseAlgorithm",
    format = Documenter.HTML(),
    modules = [PulseAlgorithm], 
    checkdocs = :none
)

deploydocs(
    repo = "github.com/EstebanLeiva/PulseAlgorithm.jl.git",
)