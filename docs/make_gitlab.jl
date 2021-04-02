using Documenter
using Marta

include("makedocs.jl")

deploydocs(
    repo = "gitlab.com/homodyne-ct/Marta.git",
    deploy_config = Documenter.GitLab()
)
