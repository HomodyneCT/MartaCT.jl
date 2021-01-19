module Interpolation

Base.Experimental.@optlevel 3

export AbstractInterpolation, AbstractBilinearInterpolation
export BilinearInterpolation
export interpolate


abstract type AbstractBilinearInterpolation end

struct BilinearInterpolation{M <: AbstractMatrix{<:Number}} <: AbstractBilinearInterpolation
    mat::M
end


@inline function (interp::BilinearInterpolation)(y::T, x::T) where {T <: Number}
    x1 = floor(Int, x)
    y1 = floor(Int, y)
    x2 = ceil(Int, x)
    y2 = ceil(Int, y)
    blerp(interp.mat, y1, y2, x1, x2, y, x)
end


@inline function interpolate(mat::AbstractMatrix{T}) where {T <: Number}
    BilinearInterpolation(mat)
end


@inline function lerp(q1::Q, q2::Q, t::T)::T where {Q <: Number,T <: Number}
    (one(T) - t) * q1 + t * q2
end

@inline function lerp(f::Function, x1::X, x2::X, x::T)::T where {X <: Number,T <: Number}
    x1 == x2 && return f(x1)
    lerp(f(x1), f(x2), (x - x1) / (x2 - x1))
end

@inline function blerp(q11::Q, q12::Q, q21::Q, q22::Q, t1::T, t2::T)::T where {Q <: Number,T <: Number}
    lerp(lerp(q11, q21, t1), lerp(q12, q22, t1), t2)
end

@inline function blerp(mat::AbstractMatrix{T}, q1::Q, q2::Q, p1::P, p2::P, q::T, p::T)::T where {Q <: Number,P <: Number,T <: Number}
    if p1 == p2
        q1 == q2 && return @inbounds mat[q1, p1]
        @inbounds q11 = mat[q1, p1]
        @inbounds q21 = mat[q2, p1]
        return lerp(q11, q21, (q - q1) / (q2 - q1))
    end
    if q1 == q2
        @inbounds q11 = mat[q1, p1]
        @inbounds q12 = mat[q1, p2]
        return lerp(q11, q12, (p - p1) / (p2 - p1))
    end
    t1 = (q - q1) / (q2 - q1)
    t2 = (p - p1) / (p2 - p1)
    @inbounds q11 = mat[q1, p1]
    @inbounds q12 = mat[q1, p2]
    @inbounds q21 = mat[q2, p1]
    @inbounds q22 = mat[q2, p2]
    blerp(q11, q12, q21, q22, t1, t2)
end

end # module
