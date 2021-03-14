module RadonAlgorithm

Base.Experimental.@optlevel 3

export Radon, RadonSquare
export AbstractFBP
export FBP, FBPFFTSquare, FBPAFFT, FBPAFFTSquare

using ..Applicative
using ..Monads
using ..CTImages
using ..Geometry
using ..FanBeam: fan2para, para2fan
import ..Marta: datatype, linspace
import ..AbstractAlgorithms:
    radon, iradon, project_image, reconstruct_image, @_alg_progress
using ..AbstractAlgorithms
using ..Interpolation: AbstractInterp2DOrNone, interpolate
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


macro defiradonalgfn(A::Symbol, f::Symbol)
    quote
        $(esc(A))(sinog::AbstractMatrix; kwargs...) =
            $f(sinog; $(esc(A)).filter, kwargs...)
        function $(esc(A))(
            sinog::AbstractMatrix,
            xs::AbstractVector,
            ys::AbstractVector,
            ϕs::ClosedInterval;
            kwargs...
        )
            $f(image, xs, ys, ϕs; $(esc(A)).filter, kwargs...)
        end
    end
end


macro defiradonfngeom(f::Symbol)
    quote
        @inline function $f(
            sinog::AbstractMatrix,
            g::AbstractParallelBeamGeometry;
            kwargs...
        )
            rows = geometry.rows
            cols = geometry.cols
            α = geometry.α
            α₀ = geometry.α₀
            $f(sinog; rows, cols, α, α₀, kwargs...)
        end
    end
end


macro radonprogress(n, p, dt=0.2)
    :(@_alg_progress "Computing Radon transform..." $n $p $dt)
end

macro iradonprogress(n, p, dt=0.2)
    :(@_alg_progress "Computing inverse Radon transform..." $n $p $dt)
end


function _make_ϕs(::Type{T}, ϕ₀, Δϕ, nϕ) where T
    range(T(ϕ₀); step = T(Δϕ / nϕ), length = nϕ)
end

function _make_ϕs(::Type{T}, ϕs::ClosedInterval, nϕ) where T
    _make_ϕs(T, minimum(ϕs), width(ϕs), nϕ)
end


abstract type AbstractFBP <: AbstractIRadonAlgorithm end

include("Filters.jl")
include("radon/radon.jl")
include("iradon/iradon.jl")

end # module
