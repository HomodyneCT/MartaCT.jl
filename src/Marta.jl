module HomodyneImaging

const _version = v"0.1.0"
version() = _version

using Reexport

include("basic_definitions.jl")
include("Applicative.jl")
include("Monads.jl")
include("TypeDict.jl")
include("CTIO.jl")
include("Info.jl")
include("Interpolation.jl")
include("Geometry.jl")
include("CTImages.jl")
include("CTData.jl")
include("FanBeam.jl")
include("Filters.jl")
include("AbstractAlgorithms.jl")
include("RadonAlgorithm.jl")
include("CalibrationBase.jl")
include("TestImages.jl")
include("Calibration.jl")
include("CTScan.jl")
include("CTPlots.jl")

@reexport using .Applicative, .Monads, .Info, .Geometry, .CTIO, .CTImages
@reexport using .AbstractAlgorithms, .RadonAlgorithm,
@reexport using .TestImages, .Calibration, .CTScan, .CTPlots

end # module
