struct Radon <: AbstractProjectionAlgorithm end

@defradonalgfn Radon radon_default


"""radon_default(image::AbstractMatrix; <keyword arguments>)

Compute the Radon transform of `image` inside a circle of
radius `hypot(rows,cols)/2` where `rows` and `cols` are the
dimensions of `image`.

See Also: [`radon_square`](@ref)
"""
function radon_default(
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
    nd = length(ts)
    nϕ = length(ϕs)
    t₀::T = hypot(rows, cols) / 2
    sθ, cθ = atan(rows, cols) |> sincos
    x₀::T = t₀ * cθ + 1
    y₀::T = t₀ * sθ + 1
    scϕs = sincos.(ϕs)
    z::T = maybe(zero(T), background)
    rmat = similar(image, nd, nϕ)
    fill!(rmat, z)
    indices = Vector{NTuple{2,Int}}(undef, length(rmat))
    tϕs = similar(indices, NTuple{2,T})
    foreach(eachindex(rmat)) do k
        ix = (k - 1) % nd
        iϕ = (k - 1) ÷ nd
        @inbounds t = ts[ix + 1]
        @inbounds s, c = scϕs[iϕ + 1]
        @inbounds indices[k] = k, iϕ
        @inbounds tϕs[k] = t * c, t * s
    end
    rimage = rescaled ? rescale(image) : image
    interp = isnothing(interpolation) ?
        interpolate(rimage) : interpolation(rimage)
    p = @radonprogress length(rmat) progress
    Threads.@threads for (k, iϕ) ∈ indices
        @inbounds tx, ty = tϕs[k]
        prex, prey = tx + x₀, ty + y₀
        o = iϕ * nd
        @inbounds rmat[k] = sum(view(tϕs, o+1:o+nd)) do (ty, tx)
            x, y = prex - tx, prey + ty
            return x ∈ 1..cols && y ∈ 1..rows ? interp(y, x) : z
        end
        next!(p)
    end
    CTSinogram(rmat)
end


@inline function radon_default(
    image::AbstractMatrix{T},
    nd::Optional{I1} = nothing,
    nϕ::Optional{I2} = nothing;
    α::Real = 360,
    α₀::Real = 0,
    ν::Real = 1,
    kwargs...
) where {
    T <: Real,
    I1 <: Integer,
    I2 <: Integer,
    U <: Real,
    Interp <: AbstractInterp2DOrNone,
}
    rows, cols = size(image)
    h::T = hypot(rows, cols)
    nd = isnothing(nd) ? round(Int, h) : nd
    nϕ = isnothing(nϕ) ? 2 * ((rows + 1) * (cols + 1) ÷ (2nd)) + 1 : nϕ
    t₀ = h / 2
    ts = linspace(-t₀..t₀, nd) * ν
    ϕ₀::T = deg2rad(α₀)
    Δϕ::T = deg2rad(α) / nϕ
    scϕs = range(ϕ₀; step = Δϕ, length = nϕ)
    radon_default(image, ts, ϕs; kwargs...)
end
