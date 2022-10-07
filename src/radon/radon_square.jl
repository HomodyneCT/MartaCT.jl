struct RadonSquare <: AbstractProjectionAlgorithm end

_alg_method(::RadonSquare) = IsRadonSquare()
@_defradonalgfn RadonSquare radon_square


"""
    radon_square(image::AbstractMatrix, ts::AbstractVector, ϕs::AbstractVector; <keyword arguments>)

Compute the Radon transform of `image` inside the circle
contained in the square of side `min(rows,cols)` where `rows` and `cols` are the
dimensions of `image`.

See also: [`radon_diag`](@ref)
"""
function radon_square end

@_defradonfn radon_square begin
    @assert 0 ∈ first(ts)..last(ts)
    x₀::T = (cols + 1) / 2
    y₀::T = (rows + 1) / 2
    κx::T = ν * (x₀ - one(Tₑ)) / half(ts) * ifelse(τ < oneunit(τ), τ, one(τ))
    κy::T = ν * (y₀ - one(Tₑ)) / half(ts) * ifelse(τ < oneunit(τ), one(τ), inv(τ))
    txs = @. T(ts * κx)
    tys = @. T(ts * κy)
    p = _radon_progress(length(scϕs), progress)
    Threads.@threads for iϕ ∈ eachindex(scϕs)
        @inbounds s, c = scϕs[iϕ]
        @inbounds @simd for j ∈ eachindex(ts)
            tx, ty = txs[j], tys[j]
            prex = tx * c + x₀
            prey = ty * s + y₀
            for k ∈ eachindex(ts)
                wx, wy = txs[k], tys[k]
                x = prex - wx * s
                y = prey + wy * c
                if 1 <= x <= cols && 1 <= y <= rows
                    sinog[j,iϕ] += interp(y, x)
                end
            end
        end
        next!(p)
    end
    δt::T = ν * width(ts) / (length(ts) - 1)
    sinog .*= δt^2
end


# function radon_square(
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
#     l = min(rows, cols)
#     nd = length(ts)
#     nt = min(l, nd)
#     nϕ = length(ϕs)
#     x′₀ = (nd - nt) ÷ 2
#     x₀::T = (cols + 1) / 2
#     y₀::T = (rows + 1) / 2
#     scϕs = sincos.(ϕs)
#     rimage = rescaled ? rescale(image) : image
#     interp = isnothing(interpolation) ?
#         interpolate(rimage) : interpolation(rimage)
#     z::T = maybe(zero(T), background)
#     rmat = similar(image, nd, nϕ)
#     fill!(rmat, z)
#     p = _radon_progress(nϕ, progress)
#     Threads.@threads for iϕ ∈ 1:nϕ
#         @inbounds s, c = scϕs[iϕ]
#         @inbounds @simd for it ∈ 1:nt
#             j = it + x′₀
#             t = ts[j]
#             prex = t * c + x₀
#             prey = t * s + y₀
#             for z ∈ ts
#                 x = prex - z * s
#                 y = prey + z * c
#                 if x ∈ 1..cols && y ∈ 1..rows
#                     rmat[j,iϕ] += interp(y,x)
#                 end
#             end
#         end
#         next!(p)
#     end
#     CTSinogram(rmat)
# end
