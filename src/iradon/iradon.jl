import .Filters: AbstractCTFilter, CTFilter, RamLak

include("fbp_default.jl")
include("fbp_square.jl")


function iradon(
    image::AbstractMatrix,
    alg::AbstractIRadonAlgorithm = FBP();
    kwargs...
)
    alg(image; kwargs...)
end

function iradon(
    image::AbstractMatrix,
    xs::AbstractVector,
    ys::AbstractVector,
    alg::AbstractIRadonAlgorithm = FBP();
    kwargs...
)
    alg(image, xs, ys; kwargs...)
end

function iradon(
    image::AbstractMatrix,
    geometry::AbstractParallelBeamGeometry
    alg::AbstractIRadonAlgorithm = FBP();
    kwargs...
)
    rows = geometry.rows
    cols = geometry.cols
    α = geometry.α
    α₀ = geometry.α₀
    radon(image, alg; nd, nϕ, α, α₀, kwargs...)
end

function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractFanBeamGeometry,
    alg::AbstractIRadonAlgorithm = FBP();
    kwargs...
)
    g′, fan_sinog = para2fan(sinog, geometry)
    iradon(fan_sinog, g′, alg; kwargs...)
end

function iradon(
    g::AbstractGeometry,
    alg::AbstractIRadonAlgorithm = FBP();
    kwargs...
)
    x -> iradon(x, g, alg; kwargs...)
end
