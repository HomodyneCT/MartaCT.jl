struct IsRadonDiag end
struct IsRadonSquare end


include("radon_default.jl")
include("radon_square.jl")


@inline function _radon(
    ::IsRadonDiag,
    a::A,
    image::AbstractMatrix{T};
    nd::Optional{I} = nothing,
    nϕ::Optional{J} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    ν::Real = 1,
    kwargs...
) where {
    A <: AbstractProjectionAlgorithm,
    T <: Real,
    I <: Integer,
    J <: Integer,
}
    rows, cols = size(image)
    h::T = hypot(rows, cols)
    nd = isnothing(nd) ? round(Int, h) : nd
    nϕ = isnothing(nϕ) ? 2 * ((rows + 1) * (cols + 1) ÷ (2nd)) + 1 : nϕ
    t₀ = h / 2
    ts = linspace(-t₀..t₀, nd) * ν
    ϕ₀::T = deg2rad(α₀)
    ϕ₁::T = ϕ₀ + deg2rad(α)
    ϕs = linspace(ORI(ϕ₀..ϕ₁), nϕ)
    a(image, ts, ϕs; kwargs...)
end


@inline function _radon(
    ::IsRadonSquare,
    a::A,
    image::AbstractMatrix{T};
    nd::Optional{I} = nothing,
    nϕ::Optional{J} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    ν::Real = 1,
    kwargs...
) where {
    A <: AbstractProjectionAlgorithm,
    T <: Real,
    I <: Integer,
    J <: Integer,
}
    rows, cols = size(image)
    l = min(rows, cols)
    nd = maybe(l, nd)
    nϕ = isnothing(nϕ) ? 2 * ((rows + 1) * (cols + 1) ÷ (2nd)) + 1 : nϕ
    t₀::T = (l - 1) / 2 * ν
    ts = linspace(-t₀..t₀, nd)
    ϕ₀::T = deg2rad(α₀)
    ϕ₁::T = ϕ₀ + deg2rad(α)
    ϕs = linspace(ORI(ϕ₀..ϕ₁), nϕ)
    a(image, ts, ϕs; kwargs...)
end


@inline function radon(
    image::AbstractMatrix,
    alg::AbstractProjectionAlgorithm = Radon();
    kwargs...
)
    alg(image; kwargs...)
end

@inline function radon(
    image::AbstractMatrix,
    ts::AbstractVector,
    ϕs::AbstractVector,
    alg::AbstractProjectionAlgorithm = Radon();
    kwargs...
)
    alg(image, ts, ϕs; kwargs...)
end

@inline function radon(
    image::AbstractMatrix,
    geometry::AbstractParallelBeamGeometry,
    alg::AbstractProjectionAlgorithm = Radon();
    kwargs...
)
    nd = geometry.nd
    nϕ = geometry.nϕ
    α = geometry.α
    α₀ = geometry.α₀
    radon(image, alg; nd, nϕ, α, α₀, kwargs...)
end

@inline function radon(
    image::AbstractMatrix,
    geometry::AbstractFanBeamGeometry,
    alg::AbstractProjectionAlgorithm = Radon();
    kwargs...
)
    g′ = ParallelBeamGeometry(geometry)
    sinog = radon(image, g′, alg; kwargs...)
    _, fan_sinog = para2fan(sinog, geometry)
    fan_sinog
end

@inline function radon(
    g::AbstractGeometry,
    alg::AbstractProjectionAlgorithm = Radon();
    kwargs...
)
    x -> radon(x, g, alg; kwargs...)
end
