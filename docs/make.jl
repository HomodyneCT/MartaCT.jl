include("makedocs.jl")

deploydocs(
    repo = "github.com/HomodyneCT/Marta.git",
    deploy_config = Documenter.GitHubActions()
)
