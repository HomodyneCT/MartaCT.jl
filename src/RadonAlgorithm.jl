module RadonAlgorithm

Base.Experimental.@optlevel 3

export Radon, RadonSquare
export FBP

using ..Applicative
using ..Monads
using ..CTImages
using ..Geometry
using ..FanBeam: fan2para, para2fan
import ..Marta: datatype, linspace
import ..AbstractAlgorithms:
    radon, iradon, project_image, reconstruct_image, @_alg_progress
using ..AbstractAlgorithms
using ..Interpolation: AbstractInterp2DOrNone, BilinearInterpolation
using FFTW, ProgressMeter, IntervalSets


macro defradonalgfn(A::Symbol, f::Symbol)
    quote
        $(esc(A))(image::AbstractMatrix; kwargs...) = $f(image; kwargs...)
        function $(esc(A))(
            image::AbstractMatrix,
            ts::AbstractVector,
            ϕs::AbstractVector;
            kwargs...
        )
            $f(image, ts, ϕs; kwargs...)
        end
    end
end


macro radonprogress(n, p, dt=0.2)
    :(@_alg_progress "Computing Radon transform..." $n $p $dt)
end

macro iradonprogress(n, p, dt=0.2)
    :(@_alg_progress "Computing inverse Radon transform..." $n $p $dt)
end


include("Filters.jl")
include("radon/radon.jl")
include("iradon/iradon.jl")

end # module
