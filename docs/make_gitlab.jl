using Documenter
using MartaCT

include("makedocs.jl")

deploydocs(
    repo = "gitlab.com/homodyne-ct/MartaCT.jl.git",
    deploy_config = Documenter.GitLab()
)
