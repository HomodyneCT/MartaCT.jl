import ..Filters: AbstractCTFilter, RamLak

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


"""
    iradon(sinog::AbstractMatrix[, algorithm::AbstractIRadonAlgorithm]; <keyword arguments>)

Compute the inverse Radon transform of `sinog` with parameters given as keyword
arguments. The parameters depend on the algorithm, please see the relative
documentation. The default algorithm is [`FBP`](@ref).

See Also: [`reconstruct_image`](@ref)
"""
@inline function iradon(
    sinog::AbstractMatrix,
    alg::AbstractIRadonAlgorithm;
    kwargs...
)
    alg(sinog; kwargs...)
end


"""
    iradon(sinog::AbstractMatrix, xs[, ys[, algorithm::AbstractIRadonAlgorithm[, coo::AbstractCoordinates]]]; <keyword arguments>)

Compute the inverse Radon transform of `sinog` on the points given by the
vectors `xs` and `ys`. If `ys` is omitted, the reconstruction is performed on
the square with `xs == ys`. The default reconstruction algorithm is
[`FBP`](@ref). Please refer to the respective documentation for additional
parameters.

See Also: [`reconstruct_image`](@ref)
"""
@inline function iradon(
    sinog::AbstractMatrix,
    xs::Union{AbstractVector,Interval},
    ys::Union{AbstractVector,Interval},
    alg::AbstractIRadonAlgorithm,
    coo::AbstractCoordinates = Cartesian();
    kwargs...
)
    alg(sinog, xs, ys, coo; kwargs...)
end


@inline function iradon(
    sinog::AbstractMatrix,
    xs::AbstractVector,
    alg::AbstractIRadonAlgorithm,
    coo::AbstractCoordinates = Cartesian();
    kwargs...
)
    alg(sinog, xs, coo; kwargs...)
end


"""
    iradon(sinog::AbstractMatrix[, geometry::AbstractGeometry[, params::AbstractParams[, algorithm::AbstractIRadonAlgorithm]]]; <keyword arguments>)

Compute the inverse Radon transform of `sinog` with explicit geometry. If the
algorithm needs specific parameters, these can be passed with `params`.
Additional parameters can be passed to the algorithm through keyword arguments.
Please see the respective documentation for more details. The default algorithm
is [`FBP`](@ref).

See Also: [`reconstruct_image`](@ref)
"""
@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractParallelBeamGeometry,
    alg::AbstractIRadonAlgorithm;
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
    alg::AbstractIRadonAlgorithm;
    kwargs...
)
    g′, para_sinog = fan2para(sinog, geometry)
    iradon(para_sinog, g′, alg; kwargs...)
end


@inline function iradon(
    g::AbstractGeometry,
    alg::AbstractIRadonAlgorithm;
    kwargs...
)
    x -> iradon(x, g, alg; kwargs...)
end


@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractParallelBeamGeometry,
    params::AbstractParams,
    alg::AbstractIRadonAlgorithm;
    kwargs...
)
    rows = geometry.rows
    cols = geometry.cols
    α = geometry.α
    α₀ = geometry.α₀
    alg(sinog, params; rows, cols, α, α₀, kwargs...)
end


@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractFanBeamGeometry,
    params::AbstractParams,
    alg::AbstractIRadonAlgorithm;
    kwargs...
)
    g′, para_sinog = fan2para(sinog, geometry)
    iradon(sinog, g′, params, alg; kwargs...)
end
