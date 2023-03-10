module RadonAlgorithm

Base.Experimental.@optlevel 3

export Radon, RadonSquare
export AbstractFBP
export FBP, FBPFFTSquare, FBPAFFT, FBPAFFTSquare
export RadonInfo, FBPInfo

using ..Applicative
using ..Monads
using ..Geometry
using ..FanBeam: fan2para, para2fan
import ..Utils: ORI, linspace, _atype, half, width
import ..AbstractAlgorithms: radon, iradon
#import ..AbstractAlgorithms: _alg_progress
using ..AbstractAlgorithms, ..Coordinates
using ..Interpolation: AbstractInterp2DOrNone, interpolate
using AbstractFFTs, IntervalSets
import FFTW


function _radon end
function _iradon end
function _alg_method end

@inline function _radon(
    a::AbstractProjectionAlgorithm,
    img::AbstractMatrix,
    args...;
    kwargs...
)
    _radon(_alg_method(a), a, img, args...; kwargs...)
end

@inline function _iradon(
    a::AbstractIRadonAlgorithm,
    sinog::AbstractMatrix,
    args...;
    kwargs...
)
    _iradon(_alg_method(a), a, sinog, args...; kwargs...)
end


macro _defradonfn(f::Symbol, body)
    esc(quote
        @inline function $f(
            image::AbstractMatrix{T},
            ts::AbstractVector,
            ϕs::AbstractVector;
            ν::Real = 1,
            scale::Optional{Real} = nothing,
            τ::Optional{Real} = nothing,
            ratio::Optional{Real} = nothing,
            background::Optional = nothing,
            interpolation::Optional{Interp} = nothing,
            progress::Bool = false,
        ) where {T,Interp <: AbstractInterp2DOrNone}
            ν = maybe(ν, scale)
            @assert ν > 0
            rows, cols = size(image)
            τ = maybe(maybe(rows / cols, τ), ratio)
            @assert τ > 0
            nd = length(ts)
            nϕ = length(ϕs)
            scϕs = map(ϕs) do ϕ
                s, c = sincos(ϕ)
                T(s), T(c)
            end
            sinog = similar(image, nd, nϕ)
            fill!(sinog, convert(T, maybe(0.0f0, background)))
            interp = isnothing(interpolation) ?
                interpolate(image) : interpolation(image)
            $body
        end
    end)
end


macro _defiradonfn(f::Symbol, body)
    esc(quote
        @inline function $f(
            sinog::AbstractMatrix,
            xs::AbstractVector,
            ys::AbstractVector,
            ::Cartesian;
            ν::Real = 1,
            scale::Optional{Real} = nothing,
            ϕs::Optional{I} = nothing,
            angles::Optional{J} = nothing,
            background::Optional = nothing,
            filter::Optional{F} = nothing,
            interpolation::Optional{Interp} = nothing,
            progress::Bool = false,
        ) where {
            I <: Interval{:closed},
            J <: Interval{:closed},
            F <: AbstractCTFilter,
            Interp <: AbstractInterp2DOrNone,
        }
            ν = maybe(ν, scale)
            @assert ν > 0
            ϕs = maybe(ϕs, angles)
            nd, nϕ = size(sinog)
            cols = length(xs)
            rows = length(ys)
            T = real(eltype(sinog))
            filtered = apply(maybe(RamLak(), filter)) do f
                filter_freq = fft(sinog, 1) .* f(T, nd)
                bfft(filter_freq, 1) |> real
            end
            interpolation = maybe(interpolate, interpolation)
            interp = interpolation(filtered)
            t₀ = (nd + 1) / 2
            ϕs = maybe(ORI(0..2π), ϕs)
            scϕs = sincos.(linspace(T, ϕs, nϕ))
            z::T = maybe(zero(T), background)
            tomog = similar(sinog, T, rows, cols)
            fill!(tomog, z)
            convert(_atype(sinog), $body)
        end
        @inline function $f(
            sinog::AbstractMatrix,
            xs::AbstractVector,
            ys::AbstractVector;
            kwargs...
        )
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
        @inline function (a::$A)(
            image::AbstractMatrix,
            ts::Interval,
            ϕs::Interval;
            kwargs...
        )
            _radon(a, image, ts, ϕs; kwargs...)
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
            xs::Interval,
            ys::Interval,
            coo::AbstractCoordinates = Cartesian();
            kwargs...
        )
            _iradon(a, sinog, xs, ys, coo; kwargs...)
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


# function _radon_progress(n::Integer, p::Bool, dt::Real=0.2)
#     _alg_progress(Progress, "Computing Radon transform...", n, p, dt)
# end

# function _iradon_progress(n::Integer, p::Bool, dt::Real=0.2)
#     _alg_progress(Progress, "Computing inverse Radon transform...", n, p, dt)
# end


abstract type AbstractFBP <: AbstractIRadonAlgorithm end

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
