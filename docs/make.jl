using Documenter
using Marta

makedocs(
    sitename = "Marta",
    format = Documenter.HTML(),
    modules = [Marta]
)

deploydocs(
    repo = "github.com/HomodyneCT/Marta.jl.git",
    deploy_config = GitHubActions()
)

deploydocs(
    repo = "gitlab.com/homodyne-ct/Marta.git",
    deploy_config = GitLab()
)
