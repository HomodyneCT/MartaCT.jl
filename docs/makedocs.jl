using Documenter
using Marta

makedocs(
    sitename = "Marta",
    format = Documenter.HTML(
        prettyurls = "local" ∉ ARGS,
    ),
    modules = [Marta],
    pages = [
        "Home" => "index.md",
        "Getting Started" => [
            "Installation.md",
        ],
        "Manual" => [
            "Tutorial.md",
        ],
        "API Reference" => [
            "api/Algorithms.md",
            "api/CTImages.md",
            "api/Geometry.md",
            "api/CTIO.md",
            "api/CTPlots.md",
            "api/Simulations.md",
            "api/Calibration.md",
            "api/TestImages.md",
            "api/CTScan.md",
            "api/Info.md",
            "api/Interpolation.md",
            "api/Applicative.md",
            "api/Monads.md",
            "api/Utils.md",
        ],
        "Index" => "GeneralIndex.md",
    ]
)
