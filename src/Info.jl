module Info

export AbstractCTInfo, CTInfo

import ..Utils: yaml_repr


abstract type AbstractCTInfo end


struct CTInfo{NT<:NamedTuple} <: AbstractCTInfo
    value::NT
end


CTInfo(pairs::Pair...) = CTInfo((; pairs...))
CTInfo(; kwargs...) = CTInfo((; kwargs...))


Base.getproperty(ctinfo::CTInfo, p::Symbol) =
    p == :value ? getfield(ctinfo, :value) :
    getproperty(getfield(ctinfo, :value), p)

Base.eltype(::Type{CTInfo{NT}}) where {NT} = eltype(NT)

Base.length(ctinfo::CTInfo) = length(ctinfo.value)
Base.iterate(ctinfo::CTInfo) = iterate(ctinfo.value)
Base.iterate(ctinfo::CTInfo, state) = iterate(ctinfo.value, state)

Base.keys(ctinfo::CTInfo) = keys(ctinfo.value)
Base.values(ctinfo::CTInfo) = values(ctinfo.value)

function Base.merge(a::CTInfo, others::CTInfo...)
    CTInfo(merge(a.value, getproperty.(others, :value)...))
end

yaml_repr(ctinfo::CTInfo) = yaml_repr(ctinfo.value)

end # module
