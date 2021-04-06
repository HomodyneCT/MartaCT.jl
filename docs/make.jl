include("makedocs.jl")

deploydocs(
    repo = "github.com/HomodyneCT/MartaCT.git",
    deploy_config = Documenter.GitHubActions()
)
