module Interpolation

Base.Experimental.@optlevel 3

export AbstractInterpolation
export AbstractInterpolation1D, AbstractInterpolation2D
export AbstractLinearInterpolation, AbstractBilinearInterpolation
export NoInterpolation
export LinearInterpolation, BilinearInterpolation
export AbstractInterp1DOrNone, AbstractInterp2DOrNone
export interpolate

using Base: @propagate_inbounds

abstract type AbstractInterpolation end
abstract type AbstractInterpolation1D <: AbstractInterpolation end
abstract type AbstractInterpolation2D <: AbstractInterpolation end
abstract type AbstractLinearInterpolation <: AbstractInterpolation1D end
abstract type AbstractBilinearInterpolation <: AbstractInterpolation2D end


struct NoInterpolation <: AbstractInterpolation end
struct BilinearInterpolation <: AbstractBilinearInterpolation end
struct LinearInterpolation <: AbstractLinearInterpolation end

const AbstractInterp1DOrNone = Union{AbstractInterpolation1D,NoInterpolation}
const AbstractInterp2DOrNone = Union{AbstractInterpolation2D,NoInterpolation}

const _interpolation_types = (
    :NoInterpolation,
    :LinearInterpolation,
    :BilinearInterpolation,
)


@propagate_inbounds function (::NoInterpolation)(a::AbstractArray, x, xs...)
    xs = (x, xs...)
    inds = round.(Int, xs)
    @boundscheck checkbounds(a, inds...)
    @inbounds a[inds...]
end

@propagate_inbounds function (::BilinearInterpolation)(mat::AbstractMatrix, y, x)

    x1 = floor(Int, AbstractFloat(x))
    y1 = floor(Int, AbstractFloat(y))
    x2 = ceil(Int, AbstractFloat(x))
    y2 = ceil(Int, AbstractFloat(y))
    blerp(mat, y1, y2, x1, x2, y, x)
end

@propagate_inbounds function (::LinearInterpolation)(v::AbstractVector, x)
    x1 = floor(Int, AbstractFloat(x))
    x2 = ceil(Int, AbstractFloat(x))
    lerp(v, x1, x2, x)
end


@propagate_inbounds function interpolate(a::AbstractArray, interp::AbstractInterpolation)
    @assert(
        a isa AbstractVecOrMat || interp isa NoInterpolation,
        "Interpolation of arrays with dimensions higher than 2 " *
        "is not currently supported"
    )
    @propagate_inbounds function (xs...)
        interp(a, xs...)
    end
end


@inline interpolate(mat::AbstractMatrix{T}) where {T <: Number} =
    interpolate(mat,  BilinearInterpolation())

@inline interpolate(v::AbstractVector{T}) where {T <: Number} =
    interpolate(v, LinearInterpolation())

for nm ∈ _interpolation_types
    @eval begin
        @inline function (interp::$nm)(a::AbstractArray{T}) where {T <: Number}
            interpolate(a, interp)
        end
    end
end


@inline function lerp(q1::Q, q2::Q, t) where {Q}
    t′::Q = t
    (oneunit(Q) - t′) * q1 + t′ * q2
end

@inline function lerp(f::Function, x1::X, x2::X, x) where {X}
    x1 == x2 && return f(x1)
    lerp(f(x1), f(x2), (x - x1) / (x2 - x1))
end

@propagate_inbounds function lerp(v::AbstractVector, x1::X, x2::X, x) where {X}
    @boundscheck checkbounds(v, x1)
    x1 == x2 && return @inbounds v[x1]
    @inbounds lerp(v[x1], v[x2], x)
end

@inline function blerp(q11::Q, q12::Q, q21::Q, q22::Q, t₁, t₂) where {Q}
    t̄₁::Q = t₁
    t̄₂::Q = t₂
    lerp(lerp(q11, q21, t̄₁), lerp(q12, q22, t̄₁), t̄₂)
end


@propagate_inbounds function blerp(
    mat::AbstractMatrix, q1::Q, q2::Q, p1::P, p2::P, q, p
) where {Q <: Integer,P <: Integer}
    @boundscheck checkbounds(mat, q1, p1)
    @boundscheck checkbounds(mat, q1, p2)
    @boundscheck checkbounds(mat, q2, p1)
    @boundscheck checkbounds(mat, q2, p2)
    @inbounds if p1 == p2
        q1 == q2 && return mat[q1, p1]
        q11 = mat[q1, p1]
        q21 = mat[q2, p1]
        lerp(q11, q21, (q - q1) / (q2 - q1))
    elseif q1 == q2
        q11 = mat[q1, p1]
        q12 = mat[q1, p2]
        lerp(q11, q12, (p - p1) / (p2 - p1))
    else
        t1 = (q - q1) / (q2 - q1)
        t2 = (p - p1) / (p2 - p1)
        q11 = mat[q1, p1]
        q12 = mat[q1, p2]
        q21 = mat[q2, p1]
        q22 = mat[q2, p2]
        blerp(q11, q12, q21, q22, t1, t2)
    end
end

end # module
