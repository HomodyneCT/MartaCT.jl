module Utils

Base.Experimental.@optlevel 3

export linspace, half
export yaml_repr, struct2dict
export ORI

import IntervalSets: width
using IntervalSets

const ORI{T} = Interval{:closed,:open,T}
ORI(i::Interval{L,R,T}) where {L,R,T} = ORI{T}(i)

function _atype end
function half end
function width end

width(xs::AbstractVector{T}) where T = width(Interval(extrema(xs)...))
width(xs::AbstractRange{T}) where T = abs(last(xs) - first(xs))
half(xs::AbstractVector{T}) where T = width(xs) / 2
half(i::Interval{L,R,T}) where {L,R,T} = width(i) / 2

_atype(::Type{T}) where {T <: AbstractArray} = T
_atype(a::AbstractArray{T}) where T = typeof(a)

@inline function linspace(start::T, stop::U, len::Integer) where {
    T <: Number,
    U <: Number,
}
    range(start, stop, length = len)
end

@inline function linspace(::Type{T}, start::X, stop::Y, len::Integer) where {
    T <: Number,
    X <: Number,
    Y <: Number,
}
    range(T(start), T(stop), length = len)
end

@inline function linspace(
    i::ClosedInterval{T},
    len::Integer,
) where {T <: Number}
    range(i, len)
end

@inline function linspace(
    ::Type{T},
    i::ClosedInterval{U},
    len::Integer,
) where {T <: Number,U <: Number}
    range(convert(ClosedInterval{T}, i), len)
end

@inline function linspace(
    i::ORI{T},
    len::Integer,
) where {T <: Number}
    range(i, len)
end

@inline function linspace(
    ::Type{T},
    i::ORI{U},
    len::Integer,
) where {T <: Number,U <: Number}
    range(convert(ORI{T}, i), len)
end

yaml_repr(obj) = obj
yaml_repr(ntup::NamedTuple) = struct2dict(ntup)

function struct2dict(obj)
    Dict(map(propertynames(obj)) do p
        field = getproperty(obj, p)
        p => yaml_repr(field)
    end)
end

@inline function _compute_radius(x::T, y::T, d::T)::T where {T <: Real}
    x == zero(T) && return abs(y) * d + one(T)
    y == zero(T) && return abs(x) * d + one(T)
    return √(x^2 + y^2) * d + one(T)
end

@inline function _compute_radius(x::T, y::T)::T where {T <: Real}
    x == 0 && return abs(y) + one(T)
    y == 0 && return abs(x) + one(T)
    return √(x^2 + y^2) + one(T)
end

@inline function _compute_angle(x::T, y::T, d::T)::T where {T <: Real}
    θ = atan(y, x)
    θ < 0 && return (θ + 2π) * d + one(T)
    θ * d + one(T)
end

@inline function _wrap_angle(θ::T, n::Integer, nh::Real)::T where {T <: Real}
    θ >= nh && return one(T)
    θ > n && return T(n)
    θ
end

end # module
