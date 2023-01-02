import .Filters: AbstractCTFilter, CTFilter, RamLak

struct IsIRadonDiag end
struct IsIRadonSquare end


include("fbp_fft.jl")
include("fbp_square.jl")
include("fbpa.jl")
include("fbpa_square.jl")


@inline function _iradon(
    ::IsIRadonDiag,
    a::A,
    sinog::AbstractMatrix;
    rows::Optional{I} = nothing,
    cols::Optional{J} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    kwargs...
) where {
    A <: AbstractIRadonAlgorithm,
    I <: Integer,
    J <: Integer,
}
    nd, nϕ = size(sinog)
    rows = maybe(cols, rows)
    rows = maybe(round(Int, nd / √2), rows)
    cols = maybe(rows, cols)
    ϕs = deg2rad(α₀)..(deg2rad(α₀) + deg2rad(α))
    sθ, cθ = sincos(atan(rows, cols))
    x₀, y₀ = cθ, sθ
    xs = linspace(real(eltype(sinog)), -x₀..x₀, cols)
    ys = linspace(real(eltype(sinog)), -y₀..y₀, rows)
    a(sinog, xs, ys; ϕs, kwargs...)
end


@inline function _iradon(
    ::IsIRadonDiag,
    a::A,
    sinog::AbstractMatrix,
    xsi::Interval,
    ysi::Interval,
    coo::Cartesian = Cartesian();
    rows::Optional{I} = nothing,
    cols::Optional{J} = nothing,
    kwargs...
) where {
    A <: AbstractIRadonAlgorithm,
    I <: Integer,
    J <: Integer,
}
    nd, nϕ = size(sinog)
    rows = maybe(cols, rows)
    rows = maybe(round(Int, nd / √2), rows)
    cols = maybe(rows, cols)
    xs = linspace(real(eltype(sinog)), xsi, cols)
    ys = linspace(real(eltype(sinog)), ysi, rows)
    a(sinog, xs, ys, coo; kwargs...)
end


@inline function _iradon(
    ::IsIRadonSquare,
    a::A,
    sinog::AbstractMatrix;
    rows::Optional{I} = nothing,
    cols::Optional{J} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    kwargs...
) where {
    A <: AbstractIRadonAlgorithm,
    I <: Integer,
    J <: Integer,
}
    nd, nϕ = size(sinog)
    rows = maybe(cols, rows)
    rows = maybe(round(Int, nd), rows)
    cols = maybe(rows, cols)
    ϕs = deg2rad(α₀)..(deg2rad(α₀) + deg2rad(α))
    sθ, cθ = sincos(atan(rows, cols))
    x₀, y₀ = cθ, sθ
    xs = linspace(real(eltype(sinog)), -x₀..x₀, cols)
    ys = linspace(real(eltype(sinog)), -y₀..y₀, rows)
    a(sinog, xs, ys; ϕs, kwargs...)
end


@inline function _iradon(
    ::IsIRadonSquare,
    a::A,
    sinog::AbstractMatrix,
    xsi::Interval,
    ysi::Interval,
    coo::Cartesian = Cartesian();
    rows::Optional{I} = nothing,
    cols::Optional{J} = nothing,
    kwargs...
) where {
    A <: AbstractIRadonAlgorithm,
    I <: Integer,
    J <: Integer,
}
    nd, nϕ = size(sinog)
    rows = maybe(cols, rows)
    rows = maybe(round(Int, nd), rows)
    cols = maybe(rows, cols)
    xs = linspace(real(eltype(sinog)), xsi, cols)
    ys = linspace(real(eltype(sinog)), ysi, rows)
    a(sinog, xs, ys, coo; kwargs...)
end


@inline function iradon(sinog::AbstractMatrix; kwargs...)
    iradon(sinog, FBP(); kwargs...)
end


@inline function iradon(
    sinog::AbstractMatrix,
    xs::Union{AbstractVector,Interval},
    ys::Union{AbstractVector,Interval};
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
