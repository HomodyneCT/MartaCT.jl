export datatype, linspace

function datatype end
function _atype end

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
    i::IntervalSets.ClosedInterval{T},
    len::Integer,
) where {T <: Number}
    a, b = IntervalSets.endpoints(i)
    linspace(a, b, len)
end

@inline function linspace(
    ::Type{T},
    i::IntervalSets.ClosedInterval{U},
    len::Integer
) where {T <: Number,U <: Number}
    a, b = IntervalSets.endpoints(i)
    linspace(T, a, b, len)
end
