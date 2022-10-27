using Documenter
using MartaCT

makedocs(
    sitename = "MartaCT",
    format = Documenter.HTML(
        prettyurls = "local" ∉ ARGS,
    ),
    modules = [MartaCT],
    pages = [
        "Home" => "index.md",
        "Getting Started" => [
            "Installation.md",
        ],
        "Manual" => [
            "Tutorial.md",
            "Quantum Optics" => "QuantumOptics.md",
        ],
        "API Reference" => [
            "api/Algorithms.md",
            "api/CTImages.md",
            "api/Geometry.md",
            "api/Simulations.md",
            "api/Calibration.md",
            "api/CTTestImages.md",
            "api/Interpolation.md",
            "api/Applicative.md",
            "api/Monads.md",
            "api/Utils.md",
        ],
        "Index" => "GeneralIndex.md",
    ]
)
