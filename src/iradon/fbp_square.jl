struct FBPFFTSquare{F<:AbstractCTFilter} <: AbstractFBP
    filter::F
end

FBPFFTSquare() = FBPFFTSquare(RamLak())
FBPFFTSquare(f::Function) = FBPFFTSquare(CTFilter(f))

_alg_method(::FBPFFTSquare) = IsIRadonSquare()
@_defiradonalgfn FBPFFTSquare fbp_fft_square

function fbp_fft_square end

@_defiradonfn fbp_fft_square begin
    κ::T = inv(ν * min(_half(xs), _half(ys)))
    xxs = xs * κ
    yys = ys * κ
    temp_images = fill(deepcopy(tomog), 1)
    Threads.resize_nthreads!(temp_images)
    p = _iradon_progress(nϕ, progress)
    Threads.@threads for iϕ ∈ eachindex(scϕs)
        sϕ, cϕ = scϕs[iϕ]
        id = Threads.threadid()
        @inbounds img = temp_images[id]
        for ix ∈ eachindex(xs), iy ∈ eachindex(ys)
            x, y = xxs[ix], yys[iy]
            t::T = (x * cϕ + y * sϕ + 1) * t₀
            if t ∈ 1..nd
                img[iy,ix] += interp(t, T(iϕ))
            end
        end
        next!(p)
    end
    foreach(temp_images) do x
        tomog .+= x
    end
    tomog
end


# function fbp_fft_square(
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
#     t₀::T = T(nd + 1) / 2
#     l = min(rows, cols)
#     x₀, y₀ = (cols - l) ÷ 2, (rows - l) ÷ 2
#     scϕs = sincos.(linspace(T, ORI(ϕs), nϕ))
#     indices = Vector{NTuple{2,Int}}(undef, l^2)
#     xys = Vector{NTuple{2,T}}(undef, l^2)
#     @inbounds @simd for k ∈ eachindex(xys)
#         ix, iy = (k - 1) ÷ l + 1, (k - 1) % l + 1
#         indices[k] = x₀ + ix, y₀ + iy
#         xys[k] = xs[ix], ys[iy]
#     end
#     z::T = maybe(zero(T), background)
#     temp_images = fill(fill(z, rows, cols), 1)
#     Threads.resize_nthreads!(temp_images)
#     p = _iradon_progress(nϕ, progress)
#     Threads.@threads for iϕ ∈ eachindex(scϕs)
#         sϕ, cϕ = scϕs[iϕ]
#         id = Threads.threadid()
#         @inbounds img = temp_images[id]
#         @inbounds @simd for k ∈ eachindex(xys)
#             ix, iy = indices[k]
#             x, y = xys[k]
#             # To be consistent with our conventions should be '+'.
#             t::T = x * cϕ + y * sϕ + t₀
#             if t ∈ 1..nd
#                 img[iy,ix] += interp(t, iϕ)
#             end
#         end
#         next!(p)
#     end
#     tomog = similar(sinog, rows, cols)
#     CTTomogram(sum!(tomog, temp_images))
# end
