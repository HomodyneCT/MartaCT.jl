struct FBPAFFT{F<:AbstractCTFilter} <: AbstractFBP
    filter::F
end

FBPAFFT() = FBPAFFT(RamLak())
FBPAFFT(f::Function) = FBPAFFT(CTFilter(f))

_alg_method(::FBPAFFT) = IsIRadonDiag()
@_defiradonalgfn FBPAFFT fbpa_fft

function fbpa_fft end

@_defiradonfn fbpa_fft begin
    xys = Vector{NTuple{2,T}}(undef, rows * cols)
    @inbounds @simd for k ∈ eachindex(xys)
        xys[k] = xs[(k - 1) ÷ rows + 1], ys[(k - 1) % rows + 1]
    end
    p = _iradon_progress(length(tomog), progress)
    Threads.@threads for k ∈ eachindex(xys)
        x, y = xys[k]
        @inbounds tomog[k] = sum(eachindex(scϕs)) do iϕ
            sϕ, cϕ = scϕs[iϕ]
            # To be consistent with our conventions should be '+'.
            t = (x * cϕ + y * sϕ + 1) * t₀
            t ∈ 1..nd ? interp(t, iϕ) : z
        end
        next!(p)
    end
    tomog
end


# function fbpa_fft(
#     sinog::AbstractMatrix{T},
#     xs::AbstractVector{U1},
#     ys::AbstractVector{U2},
#     ϕs::ClosedInterval = 0..2π;
#     background::Optional{U3} = nothing,
#     filter::Optional{F} = nothing,
#     interpolation::Optional{Interp} = nothing,
#     progress::Bool = true,
# ) where {
#     T <: Real,
#     U1 <: Real,
#     U2 <: Real,
#     U3 <: Real,
#     F <: AbstractCTFilter,
#     Interp <: AbstractInterp2DOrNone,
# }
#     nd, nϕ = size(sinog)
#     cols = length(xs)
#     rows = length(ys)
#     filtered = apply(maybe(Filters.RamLak(), filter)) do f
#         filter_freq = fft(sinog, 1) .* f(T, nd, nϕ)
#         ifft(filter_freq, 1) |> real
#     end
#     interpolation = maybe(interpolate, interpolation)
#     interp = interpolation(filtered)
#     t₀::T = (nd + 1) / 2
#     xys = Vector{NTuple{2,T}}(undef, rows * cols)
#     @inbounds @simd for k ∈ eachindex(xys)
#         xys[k] = xs[(k - 1) ÷ rows + 1], ys[(k - 1) % rows + 1]
#     end
#     scϕs = sincos.(linspace(T, ORI(ϕs), nϕ))
#     z::T = maybe(zero(T), background)
#     tomog = similar(sinog, rows, cols)
#     fill!(tomog, z)
#     p = _iradon_progress(length(tomog), progress)
#     Threads.@threads for k ∈ eachindex(xys)
#         x, y = xys[k]
#         @inbounds tomog[k] = sum(eachindex(scϕs)) do iϕ
#             sϕ, cϕ = scϕs[iϕ]
#             # To be consistent with our conventions should be '+'.
#             t = x * cϕ + y * sϕ + t₀
#             t ∈ 1..nd ? interp(t, iϕ) : z
#         end
#         next!(p)
#     end
#     CTTomogram(tomog)
# end
