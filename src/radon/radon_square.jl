struct RadonSquare <: AbstractProjectionAlgorithm end

@defradonalgfn RadonSquare radon_square


"""radon_square(
    image::AbstractMatrix{<:Real},
    ts::AbstractVector{<:Real},
    ϕs::AbstractVector{<:Real};
    <keyword arguments>
)

Compute the Radon transform of `image` inside the circle
contained in the square of side `min(rows,cols)` where `rows` and `cols` are the
dimensions of `image`.

See also: [`radon_default`](@ref)
"""
function radon_square(
    image::AbstractMatrix{T},
    ts::AbstractVector{U1},
    ϕs::AbstractVector{U2};
    background::Optional{U3} = nothing,
    rescaled::Bool = true,
    interpolation::Optional{Interp} = nothing,
    progress::Bool = true,
) where {
    T <: Real,
    U1 <: Real,
    U2 <: Real,
    U3 <: Real,
    Interp <: AbstractInterp2DOrNone,
}
    rows, cols = size(image)
    l = min(rows, cols)
    nd = length(ts)
    nt = min(l, nd)
    nϕ = length(ϕs)
    x′₀ = (nd - nt) ÷ 2
    x₀::T = (cols + 1) / 2
    y₀::T = (rows + 1) / 2
    scϕs = sincos.(ϕs)
    rimage = rescaled ? rescale(image) : image
    interp = isnothing(interpolation) ?
        interpolate(rimage) : interpolation(rimage)
    z::T = maybe(zero(T), background)
    rmat = similar(image, nd, nϕ)
    fill!(rmat, z)
    p = @radonprogress nϕ progress
    Threads.@threads for iϕ ∈ 1:nϕ
        @inbounds s, c = scϕs[iϕ]
        @inbounds @simd for it ∈ 1:nt
            j = it + x′₀
            t = ts[j]
            prex = t * c + x₀
            prey = t * s + y₀
            for z ∈ ts
                x = prex - z * s
                y = prey + z * c
                if x ∈ 1..cols && y ∈ 1..rows
                    rmat[j,iϕ] += interp(y,x)
                end
            end
        end
        next!(p)
    end
    CTSinogram(rmat)
end


@inline function radon_square(
    image::AbstractMatrix{T};
    nd::Optional{<:Integer} = nothing,
    nϕ::Optional{<:Integer} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    ν::Real = 1,
    kwargs...
) where {T <: Real}
    rows, cols = size(image)
    l = min(rows, cols)
    nd = maybe(l, nd)
    nϕ = isnothing(nϕ) ? 2 * ((rows + 1) * (cols + 1) ÷ (2nd)) + 1 : nϕ
    t₀::T = (l - 1) / 2 * ν
    ts = linspace(-t₀..t₀, nd)
    ϕ₀::T = deg2rad(α₀)
    Δϕ::T = deg2rad(α) / nϕ
    ϕs = range(ϕ₀; step = Δϕ, length = nϕ)
    radon_square(image, ts, ϕs; kwargs...)
end
