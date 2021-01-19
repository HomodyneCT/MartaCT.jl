module Filters

export AbstractCTFilter, RamLak, CTFilter, CTFilterOrFunc
export ram_lak

using FFTW: fftfreq


abstract type AbstractCTFilter end

const CTFilterOrFunc = Union{Function,AbstractCTFilter}


struct RamLak <: AbstractCTFilter end

(f::RamLak)(args...) = ram_lak(args...)


function ram_lak(::Type{T}, nd::Int, nϕ::Int) where {T<:Real}
    repeat((T ∘ abs).(fftfreq(nd)), outer = (1, nϕ))
end

ram_lak(nd::Int, nϕ::Int) = ram_lak(Float32, nd, nϕ)


struct CTFilter{F<:Function} <: AbstractCTFilter
    filter::F
end

end # module
