struct RamLak{T<:AbstractVector} <: AbstractCTFilter
    data::T
end

RamLak(a::AbstractArray{T}, dim, fs=one(real(T))) where T = RamLak(freqs(a, dim, fs))
RamLak() = RamLak(AbstractFloat[])

(f::RamLak)(::Type{T}, nd) where {T} = freqs(complex(T), nd, one(real(T)))


function ctfilter!(M̂, M, fp::CTFilterPlan{RamLak})
    P, fs = plan(fp), freqs(fp)
    mul!(M̂, P, M)
    M̂ .*= fs
    ldiv!(M, P, M̂)
end
