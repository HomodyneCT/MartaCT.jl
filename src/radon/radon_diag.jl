struct Radon <: AbstractProjectionAlgorithm end

_alg_method(::Radon) = IsRadonDiag()
@_defradonalgfn Radon radon_diag


"""
    radon_diag(image::AbstractMatrix, ts::AbstractVector, ϕs::AbstractVector; <keyword arguments>)

Compute the Radon transform of `image` inside a circle of
radius `hypot(rows,cols)/2` where `rows` and `cols` are the
dimensions of `image`.

See Also: [`radon_square`](@ref)
"""
function radon_diag end

@_defradonfn radon_diag begin
    @assert 0 ∈ first(ts)..last(ts)
    x₀::T = (cols + 1) / 2
    y₀::T = (rows + 1) / 2
    cθ::T = half(ts) * inv(√(1+τ^2))
    sθ::T = τ * cθ
    txs = @. T(ts * (x₀ - 1) / cθ)
    tys = @. T(ts * (y₀ - 1) / sθ)
    p = _radon_progress(length(scϕs), progress)
    Threads.@threads for iϕ ∈ eachindex(scϕs)
        @inbounds s, c = scϕs[iϕ]
        @inbounds @simd for i ∈ eachindex(ts)
            prex = txs[i] * c + x₀
            prey = tys[i] * s + y₀
            for j ∈ eachindex(ts)
                x = prex - txs[j] * s
                y = prey + tys[j] * c
                if x ∈ 1..cols && y ∈ 1..rows
                    sinog[i, iϕ] += interp(y, x)
                end
            end
        end
        next!(p)
    end
    γ::T = hypot(rows, cols) / min(rows, cols)
    δt::T = ν * γ * width(ts) / (length(ts) - 1)
    sinog .*= δt^2
end


# function radon_default(
#     image::AbstractMatrix{T},
#     ts::AbstractVector{X},
#     ϕs::AbstractVector{Y};
#     background::Optional{Z} = nothing,
#     rescaled::Bool = true,
#     interpolation::Optional{Interp} = nothing,
#     progress::Bool = true,
# ) where {
#     T <: Real,
#     X <: Real,
#     Y <: Real,
#     Z <: Real,
#     Interp <: AbstractInterp2DOrNone,
# }
#     rows, cols = size(image)
#     nd = length(ts)
#     nϕ = length(ϕs)
#     t₀::T = hypot(rows, cols) / 2
#     sθ, cθ = atan(rows, cols) |> sincos
#     x₀::T = t₀ * cθ + 1
#     y₀::T = t₀ * sθ + 1
#     scϕs = sincos.(ϕs)
#     z::T = maybe(zero(T), background)
#     rmat = similar(_atype(image), nd, nϕ)
#     fill!(rmat, z)
#     indices = Vector{NTuple{2,Int}}(undef, length(rmat))
#     tϕs = similar(indices, NTuple{2,T})
#     foreach(eachindex(rmat)) do k
#         ix = (k - 1) % nd
#         iϕ = (k - 1) ÷ nd
#         @inbounds t = ts[ix + 1]
#         @inbounds s, c = scϕs[iϕ + 1]
#         @inbounds indices[k] = k, iϕ
#         @inbounds tϕs[k] = t * c, t * s
#     end
#     rimage = rescaled ? rescale(image) : image
#     interp = isnothing(interpolation) ?
#         interpolate(rimage) : interpolation(rimage)
#     p = _radon_progress(length(rmat), progress)
#     Threads.@threads for (k, iϕ) ∈ indices
#         @inbounds tx, ty = tϕs[k]
#         prex, prey = tx + x₀, ty + y₀
#         o = iϕ * nd
#         @inbounds rmat[k] = sum(view(tϕs, o+1:o+nd)) do (ty, tx)
#             x, y = prex - tx, prey + ty
#             return x ∈ 1..cols && y ∈ 1..rows ? interp(y, x) : z
#         end
#         next!(p)
#     end
#     CTSinogram(rmat)
# end
