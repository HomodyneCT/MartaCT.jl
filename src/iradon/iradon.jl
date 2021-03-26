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
    sθ, cθ = sincos(atan(rows, cols))
    x₀, y₀ = cθ / ν, sθ / ν
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
    ϕ₀ = deg2rad(α₀)
    ϕ₁ = ϕ₀ + deg2rad(α)
    sθ, cθ = sincos(atan(rows, cols))
    x₀, y₀ = cθ / ν, sθ / ν
    xs = linspace(eltype(sinog), -x₀..x₀, cols)
    ys = linspace(eltype(sinog), -y₀..y₀, rows)
    a(sinog, xs, ys, ϕ₀..ϕ₁; kwargs...)
end

@inline function iradon(sinog::AbstractMatrix; kwargs...)
    iradon(sinog, FBP(); kwargs...)
end

@inline function iradon(
    sinog::AbstractMatrix,
    xs::AbstractVector,
    ys::AbstractVector;
    kwargs...
)
    iradon(sinog, xs, ys, FBP(); kwargs...)
end

@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractGeometry;
    kwargs...
)
    iradon(sinog, geometry, FBP(); kwargs...)
end

@inline function iradon(g::AbstractGeometry; kwargs...)
    iradon(g, FBP(); kwargs...)
end
