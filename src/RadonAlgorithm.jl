module RadonAlgorithm

Base.Experimental.@optlevel 3

export AbstractIRadonAlgorithm, FBP, Radon

using ..Applicative
using ..Monads
using ..CTImages
using ..Geometry
using ..FanBeam: fan2para, para2fan
import ..Marta: datatype, linspace
import ..AbstractAlgorithms:
    radon, iradon, project_image, reconstruct_image, alg_geometry, alg_params
using ..AbstractAlgorithms
using ..Interpolation: interpolate, AbstractInterp2DOrNone
using FFTW, ProgressMeter, IntervalSets


@inline _flip(x, y) = (-y, x)
@inline _flip((x, y)) = _flip(x, y)


abstract type AbstractIRadonAlgorithm <: AbstractReconstructionAlgorithm end

include("Filters.jl")
import .Filters
include("radon/radon.jl")
include("iradon/fbp.jl")

end # module
