struct FBP{F<:AbstractCTFilter} <: AbstractFBP
    filter::F
end

FBP() = FBP(RamLak())
FBP(f::Function) = FBP(CTFilter(f))

_alg_method(::FBP) = IsIRadonDiag()
@_defiradonalgfn FBP fbp_fft

function fbp_fft end

@_defiradonfn fbp_fft begin
    xys = Vector{NTuple{2,T}}(undef, rows * cols)
    δx::T = width(xs)
    δy::T = width(ys)
    h::T = hypot(δx, δy)
    λ::T = (nd - 1) / h
    κ::T = λ / ν
    @inbounds @simd for k ∈ eachindex(xys)
        x = T(xs[(k - 1) ÷ rows + 1]) * κ
        y = T(ys[(k - 1) % rows + 1]) * κ
        xys[k] = x, y
    end
    temp_images = [similar(tomog)]
    Threads.resize_nthreads!(temp_images)
    p = _iradon_progress(nϕ, progress)
    Threads.@threads for iϕ ∈ eachindex(scϕs)
        sϕ, cϕ = scϕs[iϕ]
        @inbounds img = temp_images[Threads.threadid()]
        @inbounds @simd for k ∈ eachindex(xys)
            x, y = xys[k]
            # To be consistent with our conventions should be '+'.
            t = x * cϕ + y * sϕ + t₀
            if t ∈ 1..nd
                img[k] += interp(t, iϕ)
            end
        end
        next!(p)
    end
    foreach(temp_images) do x
        tomog .+= x
    end
    γ::T = min(rows, cols) / hypot(rows, cols) * rows / cols
    δt::T = π * γ / length(scϕs) / nd * (λ^2)
    tomog .*= δt
end


# function fbp_fft(
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
#     filtered = apply(maybe(RamLak(), filter)) do f
#         filter_freq = fft(sinog, 1) .* f(T, nd, nϕ)
#         ifft(filter_freq, 1) |> real
#     end
#     interpolation = maybe(interpolate, interpolation)
#     interp = interpolation(filtered)
#     t₀::T = (nd + 1) / 2
#     scϕs = sincos.(linspace(T, ORI(ϕs), nϕ))
#     xys = Vector{NTuple{2,T}}(undef, rows * cols)
#     @inbounds @simd for k ∈ eachindex(xys)
#         xys[k] = xs[(k - 1) ÷ rows + 1], ys[(k - 1) % rows + 1]
#     end
#     z::T = maybe(zero(T), background)
#     temp_images = fill(fill(z, rows, cols), 1)
#     Threads.resize_nthreads!(temp_images)
#     p = _iradon_progress(nϕ, progress)
#     Threads.Threads.@threads for iϕ ∈ eachindex(scϕs)
#         sϕ, cϕ = scϕs[iϕ]
#         id = Threads.threadid()
#         @inbounds img = temp_images[id]
#         @inbounds @simd for k ∈ eachindex(xys)
#             x, y = xys[k]
#             # To be consistent with our conventions should be '+'.
#             t = x * cϕ + y * sϕ + t₀
#             if t ∈ 1..nd
#                 img[k] += interp(t, iϕ)
#             end
#         end
#         next!(p)
#     end
#     tomog = similar(sinog, rows, cols)
#     CTTomogram(sum!(tomog, temp_images))
# end
