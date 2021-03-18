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


@propagate_inbounds function (::NoInterpolation)(
    a::AbstractArray{T}, x::U, xs::U...
) where {T <: Number,U <: Number}
    xs = (x, xs...)
    inds = round.(Int, xs)
    @boundscheck checkbounds(a, inds...)
    @inbounds a[inds...]
end

@propagate_inbounds function (::BilinearInterpolation)(
    mat::AbstractMatrix{T}, y::U, x::U
) where {T <: Number,U <: Number}
    x1 = floor(Int, x)
    y1 = floor(Int, y)
    x2 = ceil(Int, x)
    y2 = ceil(Int, y)
    blerp(mat, y1, y2, x1, x2, y, x)
end

@propagate_inbounds function (::LinearInterpolation)(
    v::AbstractVector{T}, x::U
) where {T <: Number,U <: Number}
    x1 = floor(Int, x)
    x2 = ceil(Int, x)
    lerp(v, x1, x2, x)
end


@propagate_inbounds function interpolate(
    a::AbstractArray{T},
    interp::AbstractInterpolation,
) where {T <: Number}
    @assert(
        a isa AbstractVecOrMat || interp isa NoInterpolation,
        "Interpolation of arrays with dimensions higher than 2 " *
        "is not currently supported"
    )
    @propagate_inbounds function (xs::U...) where {U <: Number}
        interp(a, xs...)
    end
end

interpolate(mat::AbstractMatrix) = interpolate(mat,  BilinearInterpolation())
interpolate(v::AbstractVector) = interpolate(v, LinearInterpolation())

for nm ∈ _interpolation_types
    @eval begin
        @inline function (interp::$nm)(a::AbstractArray)
            interpolate(a, interp)
        end
    end
end


@inline function lerp(q1::Q, q2::Q, t::T) where {Q <: Number,T <: Number}
    (one(Q) - t) * q1 + t * q2
end

@inline function lerp(
    f::Function, x1::X, x2::X, x::T
) where {X <: Number,T <: Number}
    x1 == x2 && return f(x1)
    lerp(f(x1), f(x2), (x - x1) / (x2 - x1))
end

@propagate_inbounds function lerp(
    v::AbstractVector{T}, x1::X, x2::X, x::U
) where {T <: Number,X <: Integer,U<:Number}
    @boundscheck checkbounds(v, x1)
    x1 == x2 && return @inbounds v[x1]
    @boundscheck checkbounds(v, x2)
    @inbounds lerp(v[x1], v[x2], x)
end

@inline function blerp(
    q11::Q, q12::Q, q21::Q, q22::Q, t1::T1, t2::T2
) where {Q <: Number,T1 <: Number, T2 <: Number}
    lerp(lerp(q11, q21, t1), lerp(q12, q22, t1), t2)
end

@propagate_inbounds function blerp(
    mat::AbstractMatrix{T}, q1::Q, q2::Q, p1::P, p2::P, q::U, p::U
) where {Q <: Integer,P <: Integer,T <: Number,U <: Number}
    @boundscheck checkbounds(mat, q1, p1)
    @boundscheck checkbounds(mat, q1, p2)
    @boundscheck checkbounds(mat, q2, p1)
    @boundscheck checkbounds(mat, q2, p2)
    @inbounds if p1 == p2
        q1 == q2 && return mat[q1, p1]
        q11 = mat[q1, p1]
        q21 = mat[q2, p1]
        return lerp(q11, q21, (q - q1) / (q2 - q1))
    end
    @inbounds if q1 == q2
        q11 = mat[q1, p1]
        q12 = mat[q1, p2]
        return lerp(q11, q12, (p - p1) / (p2 - p1))
    end
    t1 = (q - q1) / (q2 - q1)
    t2 = (p - p1) / (p2 - p1)
    @inbounds begin
        q11 = mat[q1, p1]
        q12 = mat[q1, p2]
        q21 = mat[q2, p1]
        q22 = mat[q2, p2]
        blerp(q11, q12, q21, q22, t1, t2)
    end
end

end # module
