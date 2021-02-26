module RadonAlgorithm

Base.Experimental.@optlevel 3

export AbstractIRadonAlgorithm, FBP, Radon


using ..Monads
using ..Filters
using ..CTImages
using ..Geometry
using ..FanBeam: fan2para, para2fan
import ..Marta: datatype
import ..AbstractAlgorithms:
    radon, iradon, project_image, reconstruct_image, alg_geometry, alg_params
using ..AbstractAlgorithms
using ..Interpolation: interpolate, AbstractInterp2DOrNone
using FFTW, ProgressMeter


@inline _flip(x, y) = (-y, x)
@inline _flip((x, y)) = _flip(x, y)


abstract type AbstractIRadonAlgorithm <: AbstractReconstructionAlgorithm end


include("radon/radon.jl")
include("iradon/fbp.jl")

end # module
