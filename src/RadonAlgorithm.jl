module RadonAlgorithm

Base.Experimental.@optlevel 3

export Radon, RadonSquare
export AbstractFBP
export FBP, FBPFFTSquare, FBPAFFT, FBPAFFTSquare
export RadonInfo, FBPInfo

using ..Applicative
using ..Monads
using ..CTImages
using ..Geometry
using ..FanBeam: fan2para, para2fan
import ..Utils: ORI, linspace, _atype, _half
import ..AbstractAlgorithms:
    radon,
    iradon,
    project_image,
    reconstruct_image,
    _alg_progress
using ..AbstractAlgorithms, ..Coordinates
using ..Interpolation: AbstractInterp2DOrNone, interpolate
using FFTW, ProgressMeter, IntervalSets


function _radon end
function _iradon end
function _alg_method end

@inline function _radon(
    a::AbstractProjectionAlgorithm,
    x::AbstractMatrix;
    kwargs...
)
    _radon(_alg_method(a), a, x; kwargs...)
end

@inline function _iradon(
    a::AbstractIRadonAlgorithm,
    x::AbstractMatrix;
    kwargs...
)
    _iradon(_alg_method(a), a, x; kwargs...)
end


macro _defradonfn(f::Symbol, body)
    esc(quote
        @inline function $f(
            image::AbstractMatrix{T},
            ts::AbstractVector{X},
            ϕs::AbstractVector{Y};
            ν::Real = 1,
            background::Optional{Z} = nothing,
            rescaled::Bool = false,
            interpolation::Optional{Interp} = nothing,
            progress::Bool = true,
        ) where {
            T <: Real,
            X <: Real,
            Y <: Real,
            Z <: Real,
            Interp <: AbstractInterp2DOrNone,
        }
            rows, cols = size(image)
            nd = length(ts)
            nϕ = length(ϕs)
            scϕs = sincos.(ϕs)
            z::T = maybe(zero(T), background)
            sinog = similar(_atype(image), nd, nϕ)
            fill!(sinog, z)
            rimage = rescaled ? rescale(image) : image
            interp = isnothing(interpolation) ?
                interpolate(rimage) : interpolation(rimage)
            CTSinogram($body)
        end
    end)
end


macro _defiradonfn(f::Symbol, body)
    esc(quote
        @inline function $f(
            sinog::AbstractMatrix{T},
            xs::AbstractVector{U1},
            ys::AbstractVector{U2},
            ::Cartesian;
            ν::Real = 1,
            ϕs::Optional{I} = nothing,
            background::Optional{U3} = nothing,
            filter::Optional{F} = nothing,
            interpolation::Optional{Interp} = nothing,
            progress::Bool = true,
        ) where {
            T <: Real,
            U1 <: Real,
            U2 <: Real,
            I <: Interval{:closed},
            U3 <: Real,
            F <: AbstractCTFilter,
            Interp <: AbstractInterp2DOrNone,
        }
            nd, nϕ = size(sinog)
            cols = length(xs)
            rows = length(ys)
            filtered = apply(maybe(RamLak(), filter)) do f
                filter_freq = fft(sinog, 1) .* f(T, nd, nϕ)
                ifft(filter_freq, 1) |> real
            end
            interpolation = maybe(interpolate, interpolation)
            interp = interpolation(filtered)
            t₀::T = (nd + 1) / 2
            ϕs = maybe(ORI(0..2π), ϕs)
            scϕs = sincos.(linspace(T, ϕs, nϕ))
            z::T = maybe(zero(T), background)
            tomog = similar(_atype(sinog), rows, cols)
            fill!(tomog, z)
            CTTomogram($body)
        end
        @inline function $f(
            sinog::AbstractMatrix{T},
            xs::AbstractVector{X},
            ys::AbstractVector{Y};
            kwargs...
        ) where {T<:Real,X<:Real,Y<:Real}
            $f(sinog, xs, ys, Cartesian(); kwargs...)
        end
    end)
end


macro _defradonalgfn(A::Symbol, f::Symbol)
    esc(quote
        @inline function (a::$A)(
            image::AbstractMatrix,
            ts::AbstractVector,
            ϕs::AbstractVector;
            kwargs...
        )
            $f(image, ts, ϕs; kwargs...)
        end
        @inline function (a::$A)(image::AbstractMatrix; kwargs...)
            _radon(a, image; kwargs...)
        end
    end)
end


macro _defiradonalgfn(A::Symbol, f::Symbol)
    esc(quote
        @inline function (a::$A)(
            sinog::AbstractMatrix,
            xs::AbstractVector,
            ys::AbstractVector,
            coo::AbstractCoordinates = Cartesian();
            kwargs...
        )
            $f(sinog, xs, ys, coo; kwargs...)
        end
        @inline function (a::$A)(
            sinog::AbstractMatrix,
            xs::AbstractVector,
            coo::AbstractCoordinates = Cartesian();
            kwargs...
        )
            $f(sinog, xs, xs, coo; kwargs...)
        end
        @inline function (a::$A)(sinog::AbstractMatrix; kwargs...)
            _iradon(a, sinog; kwargs...)
        end
    end)
end


function _radon_progress(n::Integer, p::Bool, dt::Real=0.2)
    _alg_progress(Progress, "Computing Radon transform...", n, p, dt)
end

function _iradon_progress(n::Integer, p::Bool, dt::Real=0.2)
    _alg_progress(Progress, "Computing inverse Radon transform...", n, p, dt)
end


abstract type AbstractFBP <: AbstractIRadonAlgorithm end

include("Filters.jl")
include("radon/radon.jl")
include("iradon/iradon.jl")


struct RadonInfo{
    G <: AbstractGeometry,
    A <: AbstractProjectionAlgorithm
} <: AlgorithmInfo{A}
    geometry::G
    algorithm::A
end

RadonInfo(g::AbstractGeometry) = RadonInfo(g, Radon())


struct FBPInfo{G<:AbstractGeometry,A<:AbstractFBP} <: AlgorithmInfo{A}
    geometry::G
    algorithm::A
end

FBPInfo(g::AbstractGeometry) = FBPInfo(g, FBP())

end # module
