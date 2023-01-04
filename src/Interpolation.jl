module Interpolation

Base.Experimental.@optlevel 3

export AbstractInterpolation
export AbstractInterpolation1D, AbstractInterpolation2D
export AbstractLinearInterpolation, AbstractBilinearInterpolation
export NoInterpolation
export LinearInterpolation, BilinearInterpolation
export AbstractInterp1DOrNone, AbstractInterp2DOrNone
export InterpolatedArray
export interpolate

using Base: @propagate_inbounds
import Base: size, getindex

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


struct InterpolatedArray{I<:AbstractInterpolation,T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
    data::A
end

function InterpolatedArray(a::AbstractArray, interp::AbstractInterpolation)
    InterpolatedArray{typeof(interp),eltype(a),ndims(a),typeof(a)}(a)
end


size(interp::InterpolatedArray) = size(interp.data)


function interpolate end
function blerp end
function lerp end


_floori(x) = unsafe_trunc(Int, floor(x))
_ceili(x) = unsafe_trunc(Int, ceil(x))
_roundi(x) = unsafe_trunc(Int, round(x))


@propagate_inbounds function getindex(iarr::InterpolatedArray{NoInterpolation}, x, xs...)
    a = iarr.data
    inds = _roundi(x), (_roundi(x) for x in xs)...
    @boundscheck checkbounds(a, inds...)
    @inbounds a[inds...]
end


@propagate_inbounds function getindex(iarr::InterpolatedArray{LinearInterpolation}, x)
    v = iarr.data
    x1 = _floori(x)
    x2 = _ceili(x)
    @boundscheck begin
        checkbounds(v, x1)
        checkbounds(v, x2)
    end
    @inbounds lerp(v, x1, x2, x)
end


@propagate_inbounds function getindex(iarr::InterpolatedArray{BilinearInterpolation,T,2}, y, x) where T
    mat = iarr.data
    q, p = y, x
    q1 = _floori(q)
    p1 = _floori(p)
    q2 = _ceili(q)
    p2 = _ceili(p)
    @boundscheck begin
        checkbounds(mat, q1, p1)
        checkbounds(mat, q1, p2)
        checkbounds(mat, q2, p1)
        checkbounds(mat, q2, p2)
    end
    @inbounds blerp(mat, q1, q2, p1, p2, q, p)
end


@propagate_inbounds function getindex(iarr::InterpolatedArray{BilinearInterpolation}, y, x, idxs::Int...)
    arr = iarr.data
    q, p = y, x
    q1 = _floori(q)
    p1 = _floori(p)
    q2 = _ceili(q)
    p2 = _ceili(p)
    @boundscheck begin
        checkbounds(arr, q1, p1, idxs...)
        checkbounds(arr, q1, p2, idxs...)
        checkbounds(arr, q2, p1, idxs...)
        checkbounds(arr, q2, p2, idxs...)
    end
    @inbounds blerp(arr, q1, q2, p1, p2, q, p, idxs...)
end


@propagate_inbounds function (iarr::InterpolatedArray{NoInterpolation})(x, xs...)
    a = iarr.data
    inds = round(Int, real(x)), (round(Int, real(x′)) for x′ in xs)...
    @boundscheck checkbounds(a, inds...)
    @inbounds a[inds...]
end


@propagate_inbounds function (iarr::InterpolatedArray{LinearInterpolation})(x)
    v = iarr.data
    x′ = real(x)
    x1 = _floori(x′)
    x2 = _ceili(x′)
    @boundscheck begin
        checkbounds(v, x1)
        checkbounds(v, x2)
    end
    @inbounds lerp(v, x1, x2, x)
end


@propagate_inbounds function (iarr::InterpolatedArray{BilinearInterpolation,T,2})(y, x) where T
    mat = iarr.data
    q, p = real(y), real(x)
    q1 = _floori(q)
    p1 = _floori(p)
    q2 = _ceili(q)
    p2 = _ceili(p)
    @boundscheck begin
        checkbounds(mat, q1, p1)
        checkbounds(mat, q1, p2)
        checkbounds(mat, q2, p1)
        checkbounds(mat, q2, p2)
    end
    @inbounds blerp(mat, q1, q2, p1, p2, q, p)
end


@propagate_inbounds function (iarr::InterpolatedArray{BilinearInterpolation})(y, x, idxs::Int...)
    mat = iarr.data
    q, p = real(y), real(x)
    q1 = _floori(q)
    p1 = _floori(p)
    q2 = _ceili(q)
    p2 = _ceili(p)
    @boundscheck begin
        checkbounds(mat, q1, p1, idxs...)
        checkbounds(mat, q1, p2, idxs...)
        checkbounds(mat, q2, p1, idxs...)
        checkbounds(mat, q2, p2, idxs...)
    end
    @inbounds blerp(mat, q1, q2, p1, p2, q, p, idxs...)
end


function interpolate(a::AbstractArray, interp::AbstractInterpolation)
    @assert(
        a isa AbstractVecOrMat || interp isa NoInterpolation || interp isa BilinearInterpolation,
        "Interpolation of arrays with dimensions higher than 2 " *
        "is not currently supported, besides for " *
        "BilinearInterpolation on the first 2 indices"
    )
    InterpolatedArray(a, interp)
end


interpolate(mat::AbstractMatrix) = interpolate(mat, BilinearInterpolation())
interpolate(v::AbstractVector) = interpolate(v, LinearInterpolation())


for nm ∈ _interpolation_types
    @eval begin
        @inline function (interp::$nm)(a::AbstractArray)
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


@propagate_inbounds function lerp(v::AbstractVector, x1::Int, x2::Int, x)
    @boundscheck begin
        checkbounds(v, x1)
        checkbounds(v, x2)
    end
    @inbounds begin
        x1 == x2 && return v[x1]
        lerp(v[x1], v[x2], x)
    end
end


@inline function blerp(q11::Q, q12::Q, q21::Q, q22::Q, t₁, t₂) where {Q}
    t̄₁::Q = t₁
    t̄₂::Q = t₂
    lerp(lerp(q11, q21, t̄₁), lerp(q12, q22, t̄₁), t̄₂)
end


@propagate_inbounds function blerp(
    mat::AbstractMatrix, q1::Int, q2::Int, p1::Int, p2::Int, q, p
)
    @boundscheck begin
        checkbounds(mat, q1, p1)
        checkbounds(mat, q1, p2)
        checkbounds(mat, q2, p1)
        checkbounds(mat, q2, p2)
    end
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


@propagate_inbounds function blerp(
    arr::AbstractArray, q1::Int, q2::Int, p1::Int, p2::Int, q, p, idxs::Int...
)
    @boundscheck begin
        checkbounds(arr, q1, p1, idxs...)
        checkbounds(arr, q1, p2, idxs...)
        checkbounds(arr, q2, p1, idxs...)
        checkbounds(arr, q2, p2, idxs...)
    end
    mat = @view arr[:,:,idxs...]
    @inbounds blerp(mat, q1, q2, p1, p2, q, p)
end

end # module
