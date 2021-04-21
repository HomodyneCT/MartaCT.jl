include("makedocs.jl")

deploydocs(
    repo = "github.com/HomodyneCT/MartaCT.jl.git",
    deploy_config = Documenter.GitHubActions()
)
