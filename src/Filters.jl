module Filters

export AbstractCTFilter, RamLak, CTFilter, CTFilterOrFunc
export ram_lak

using FFTW: fftfreq


abstract type AbstractCTFilter end

const CTFilterOrFunc = Union{Function,AbstractCTFilter}


struct RamLak <: AbstractCTFilter end

(f::RamLak)(args...) = ram_lak(args...)


function ram_lak(::Type{T}, nd::Int) where {T}
    abs.(fftfreq(nd, one(real(T))))
end

ram_lak(nd::Int) = ram_lak(Float32, nd)


struct CTFilter{F<:Function} <: AbstractCTFilter
    filter::F
end

end # module
