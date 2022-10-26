module Utils

Base.Experimental.@optlevel 3

export ORI, half, linspace

import IntervalSets: width
using IntervalSets

const ORI{T} = Interval{:closed,:open,T}
ORI(i::Interval{L,R,T}) where {L,R,T} = ORI{T}(i)

function _atype end
function half end
function linspace end

width(xs::AbstractVector{T}) where T = width(Interval(extrema(xs)...))
width(xs::AbstractRange{T}) where T = abs(last(xs) - first(xs))
half(xs::AbstractVector{T}) where T = width(xs) / 2
half(i::Interval{L,R,T}) where {L,R,T} = width(i) / 2

_atype(::Type{T}) where {T <: AbstractArray} = T
_atype(a::AbstractArray{T}) where T = typeof(a)

@inline function linspace(start, stop, len::Integer)
    range(start, stop, length = len)
end

@inline function linspace(::Type{T}, start, stop, len::Integer) where {T}
    range(T(start), T(stop), length = len)
end

@inline linspace(i::ClosedInterval, len::Integer) = range(i, len)

@inline function linspace(::Type{T}, i::ClosedInterval, len::Integer) where {T}
    range(convert(ClosedInterval{T}, i), len)
end

@inline linspace(i::ORI, len::Integer) = range(i, len)

@inline linspace(::Type{T}, i::ORI, len::Integer) where {T} = range(convert(ORI{T}, i), len)

dict_repr(obj) = obj
dict_repr(ntup::NamedTuple) = struct2dict(ntup)

function struct2dict(obj)
    Dict(map(propertynames(obj)) do p
        field = getproperty(obj, p)
        p => dict_repr(field)
    end)
end

@inline function _compute_radius(x::T, y::T, d::T)::T where {T}
    x == zero(T) && return abs(y) * d + oneunit(T)
    y == zero(T) && return abs(x) * d + oneunit(T)
    return √(x^2 + y^2) * d + oneunit(T)
end

@inline function _compute_radius(x::T, y::T)::T where {T}
    x == zero(T) && return abs(y) + one(T)
    y == zero(T) && return abs(x) + one(T)
    return √(x^2 + y^2) + one(T)
end

@inline function _compute_angle(x::T, y::T, d::T)::T where {T}
    θ = atan(y, x)
    θ < zero(T) && return (θ + 2π) * d + oneunit(T)
    θ * d + oneunit(T)
end

@inline function _wrap_angle(θ::T, n::Integer)::T where {T}
    θ >= n + 1//2 && return oneunit(T)
    θ > n && return T(n)
    θ
end

end # module
