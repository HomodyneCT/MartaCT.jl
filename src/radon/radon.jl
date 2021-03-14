include("radon_default.jl")
include("radon_square.jl")


function radon(
    image::AbstractMatrix,
    alg::AbstractProjectionAlgorithm = Radon();
    kwargs...
)
    alg(image; kwargs...)
end

function radon(
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

function radon(
    g::AbstractGeometry,
    alg::AbstractProjectionAlgorithm = Radon();
    kwargs...
)
    x -> radon(x, g, alg; kwargs...)
end
