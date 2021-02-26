export datatype, linspace

function datatype end
function _atype end

function linspace(start::T, stop::U, len::Integer) where {
    T <: Number,
    U <: Number,
}
    range(start, stop, length = len)
end

function linspace(::Type{T}, start::X, stop::Y, len::Integer) where {
    T <: Number,
    X <: Number,
    Y <: Number,
}
    range(T(start), T(stop), length = len)
end
