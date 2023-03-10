module Filters

export AbstractCTFilter, CTFilterPlan
export ctfilter!

import AbstractFFTs: Plan as FFTPlan
using AbstractFFTs
using LinearAlgebra


abstract type AbstractCTFilter end


function prepare_plan(::Type{F}, M::AbstractArray{T}, dims, fs=nothing) where {F<:AbstractCTFilter,T<:Real}
    fs′ = something(fs, one(T))
    F(M, dims, fs′), plan_rfft(M, dims)
end


function prepare_plan(::Type{F}, M::AbstractArray{T}, dims, fs=nothing) where {F<:AbstractCTFilter,T<:Complex}
    fs′ = something(fs, one(real(T)))
    F(M, dims, fs′), plan_fft(M, dims)
end


struct CTFilterPlan{F<:AbstractCTFilter,P<:FFTPlan}
    filter::F
    plan::P
end


function CTFilterPlan(::Type{F}, M::AbstractArray, dims, fs=nothing) where {F<:AbstractCTFilter}
    filter, plan = prepare_plan(F, M, dims, fs)
    CTFilterPlan(filter, plan)
end


freqs(a::AbstractArray, d::Int, fs) = freqs(eltype(a), size(a, d), fs)
freqs(::Type{T}, nd::Int, fs) where {T<:Real} = rfftfreq(nd, fs)
freqs(::Type{T}, nd::Int, fs) where {T<:Complex} = abs.(fftfreq(nd, fs))

freqs(filter::AbstractCTFilter) = filter.data
freqs(fp::CTFilterPlan) = freqs(fp.filter)

plan(fp::CTFilterPlan) = fp.plan


include("filters/ram_lak.jl")

end # module
