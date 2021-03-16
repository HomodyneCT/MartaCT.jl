module Utils

Base.Experimental.@optlevel 3


export linspace
export yaml_repr, struct2dict
export ORI

using IntervalSets


const ORI{T} = Interval{:closed,:open,T}
ORI(i::Interval) = ORI{eltype(i)}(i)


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

end # module
