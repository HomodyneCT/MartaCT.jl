import .Filters: AbstractCTFilter, CTFilter, RamLak

struct IsIRadonDiag end
struct IsIRadonSquare end


include("fbp_default.jl")
include("fbp_square.jl")
include("fbpa.jl")
include("fbpa_square.jl")


@inline function _iradon(
    ::IsIRadonDiag,
    a::A,
    sinog::AbstractMatrix{T};
    rows::Optional{I} = nothing,
    cols::Optional{J} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    ν::Real = 1,
    kwargs...
) where {
    A <: AbstractIRadonAlgorithm,
    T <: Real,
    I <: Integer,
    J <: Integer,
}
    nd, nϕ = size(sinog)
    rows = maybe(cols, rows)
    rows = maybe(round(Int, nd / √2), rows)
    cols = maybe(rows, cols)
    ϕ₀ = deg2rad(α₀)
    ϕ₁ = ϕ₀ + deg2rad(α)
    t₀ = (nd + 1) / 2
    sθ, cθ = sincos(atan(rows, cols))
    x₀, y₀ = t₀ * cθ / ν, t₀ * sθ / ν
    xs = linspace(eltype(sinog), -x₀..x₀, cols)
    ys = linspace(eltype(sinog), -y₀..y₀, rows)
    a(sinog, xs, ys, ϕ₀..ϕ₁; kwargs...)
end


@inline function _iradon(
    ::IsIRadonSquare,
    a::A,
    sinog::AbstractMatrix{T};
    rows::Optional{I} = nothing,
    cols::Optional{J} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    ν::Real = 1,
    kwargs...
) where {
    A <: AbstractIRadonAlgorithm,
    T <: Real,
    I <: Integer,
    J <: Integer,
}
    nd, nϕ = size(sinog)
    rows = maybe(cols, rows)
    rows = maybe(round(Int, nd), rows)
    cols = maybe(rows, cols)
    l = min(rows, cols)
    ϕ₀ = deg2rad(α₀)
    ϕ₁ = ϕ₀ + deg2rad(α)
    t₀ = (nd + 1) / 2ν
    xs = linspace(eltype(sinog), -t₀..t₀, l)
    ys = linspace(eltype(sinog), -t₀..t₀, l)
    a(sinog, xs, ys, ϕ₀..ϕ₁; kwargs...)
end


@inline function iradon(
    sinog::AbstractMatrix,
    alg::AbstractIRadonAlgorithm = FBP();
    kwargs...
)
    alg(sinog; kwargs...)
end

@inline function iradon(
    sinog::AbstractMatrix,
    xs::AbstractVector,
    ys::AbstractVector,
    alg::AbstractIRadonAlgorithm = FBP();
    kwargs...
)
    alg(sinog, xs, ys; kwargs...)
end

@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractParallelBeamGeometry,
    alg::AbstractIRadonAlgorithm = FBP();
    kwargs...
)
    rows = geometry.rows
    cols = geometry.cols
    α = geometry.α
    α₀ = geometry.α₀
    iradon(sinog, alg; rows, cols, α, α₀, kwargs...)
end

@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractFanBeamGeometry,
    alg::AbstractIRadonAlgorithm = FBP();
    kwargs...
)
    g′, para_sinog = fan2para(sinog, geometry)
    iradon(para_sinog, g′, alg; kwargs...)
end

@inline function iradon(
    g::AbstractGeometry,
    alg::AbstractIRadonAlgorithm = FBP();
    kwargs...
)
    x -> iradon(x, g, alg; kwargs...)
end
