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

macro _checkbounds end


macro _checkbounds(s, a, x, xs...)
    if s === :NoInterpolation
        quote
            x_ = $(esc(x))
            inds = round(Int, $x), round.(Int, $(xs...))...
            @boundscheck checkbounds($(esc(a)), inds...)
            inds
        end
    elseif s === :LinearInterpolation
        quote
            v = $(esc(a))
            x_ = $(esc(x))
            x1 = floor(Int, x_)
            x2 = ceil(Int, x_)
            @boundscheck begin
                checkbounds(v, x1)
                checkbounds(v, x2)
            end
            x1, x2
        end
    elseif s === :BilinearInterpolation
        @assert length(xs) >= 1
        quote
            c_, idxs... = $(esc.(xs)...)
            r_ = $(esc(x))
            mat = $(esc(a))
            r1 = floor(Int, r_)
            c1 = floor(Int, c_)
            r2 = ceil(Int, r_)
            c2 = ceil(Int, c_)
            @boundscheck begin
                checkbounds(mat, r1, c1, idxs...)
                checkbounds(mat, r1, c2, idxs...)
                checkbounds(mat, r2, c1, idxs...)
                checkbounds(mat, r2, c2, idxs...)
            end
            r1, r2, c1, c2
        end
    else
        error("Unknown interpolation kind, got: '$s'")
    end
end


@propagate_inbounds function getindex(iarr::InterpolatedArray{NoInterpolation}, x, xs...)
    a = iarr.data
    inds = @_checkbounds NoInterpolation a x xs...
    @inbounds a[inds...]
end


@propagate_inbounds function getindex(iarr::InterpolatedArray{LinearInterpolation}, x)
    v = iarr.data
    x1, x2 = @_checkbounds LinearInterpolation v x
    @inbounds lerp(v, x1, x2, x)
end


@propagate_inbounds function getindex(iarr::InterpolatedArray{BilinearInterpolation,T,2}, y, x) where T
    mat = iarr.data
    y1, y2, x1, x2 = @_checkbounds BilinearInterpolation mat y x
    @inbounds blerp(mat, y1, y2, x1, x2, y, x)
end


@propagate_inbounds function getindex(iarr::InterpolatedArray{BilinearInterpolation}, y, x, idxs::Int...)
    arr = iarr.data
    y1, y2, x1, x2 = @_checkbounds BilinearInterpolation arr y x idxs...
    @inbounds blerp(arr, y1, y2, x1, x2, y, x, idxs...)
end


@propagate_inbounds function (iarr::InterpolatedArray{NoInterpolation})(x, xs...)
    a = iarr.data
    inds = round(Int, x), (round(Int, x′) for x′ in xs)...
    @boundscheck checkbounds(a, inds...)
    @inbounds a[inds...]
end


@propagate_inbounds function (iarr::InterpolatedArray{LinearInterpolation})(x)
    v = iarr.data
    x1 = floor(Int, real(x))
    x2 = ceil(Int, real(x))
    @boundscheck begin
        checkbounds(v, x1)
        checkbounds(v, x2)
    end
    @inbounds lerp(v, x1, x2, x)
end


@propagate_inbounds function (iarr::InterpolatedArray{BilinearInterpolation,T,2})(y, x) where T
    mat = iarr.data
    x1 = floor(Int, real(x))
    y1 = floor(Int, real(y))
    x2 = ceil(Int, real(x))
    y2 = ceil(Int, real(y))
    @boundscheck begin
        checkbounds(mat, y1, x1)
        checkbounds(mat, y1, x2)
        checkbounds(mat, y2, x1)
        checkbounds(mat, y2, x2)
    end
    @inbounds blerp(mat, y1, y2, x1, x2, y, x)
end


@propagate_inbounds function (iarr::InterpolatedArray{BilinearInterpolation})(y, x, idxs::Int...)
    mat = iarr.data
    x1 = floor(Int, real(x))
    y1 = floor(Int, real(y))
    x2 = ceil(Int, real(x))
    y2 = ceil(Int, real(y))
    @boundscheck begin
        checkbounds(mat, y1, x1, idxs...)
        checkbounds(mat, y1, x2, idxs...)
        checkbounds(mat, y2, x1, idxs...)
        checkbounds(mat, y2, x2, idxs...)
    end
    @inbounds blerp(mat, y1, y2, x1, x2, y, x, idxs...)
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
