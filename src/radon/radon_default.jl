struct Radon <: AbstractProjectionAlgorithm end

_alg_method(::Radon) = IsRadonDiag()
@_defradonalgfn Radon radon_default


"""
    radon_default(image::AbstractMatrix, ts::AbstractVector, ϕs::AbstractVector; <keyword arguments>)

Compute the Radon transform of `image` inside a circle of
radius `hypot(rows,cols)/2` where `rows` and `cols` are the
dimensions of `image`.

See Also: [`radon_square`](@ref)
"""
function radon_default end

@_defradonfn radon_default begin
    @assert 0 ∈ first(ts)..last(ts)
    t₀::T = hypot(rows, cols) / 2
    κ::T = _half(ts) / (t₀ * ν)
    sθ, cθ = atan(rows, cols) |> sincos
    x₀::T = t₀ * cθ + 1
    y₀::T = t₀ * sθ + 1
    zs = @. T(ts / κ)
    p = _radon_progress(length(scϕs), progress)
    Threads.@threads for iϕ ∈ eachindex(scϕs)
        @inbounds s, c = scϕs[iϕ]
        @inbounds @simd for it ∈ eachindex(zs)
            t = zs[it]
            prex = t * c + x₀
            prey = t * s + y₀
            for z ∈ zs
                x = prex - z * s
                y = prey + z * c
                if x ∈ 1..cols && y ∈ 1..rows
                    sinog[it, iϕ] += interp(y, x)
                end
            end
        end
        next!(p)
    end
    sinog .*= κ^2
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
