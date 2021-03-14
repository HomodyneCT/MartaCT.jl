import .Filters: AbstractCTFilter, CTFilter, RamLak

macro defdiagkwfn(f::Symbol)
    quote
        @inline function $f(
            sinog::AbstractMatrix;
            rows::Optional{Integer} = nothing,
            cols::Optional{Integer} = nothing,
            α::Real = 360,
            α₀::Real = 0,
            ν::Real = 1,
            kwargs...
        )
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
            $f(sinog, xs, ys, ϕ₀..ϕ₁; kwargs...)
        end
    end
end

macro defsquarekwfn(f::Symbol)
    quote
        @inline function $f(
            sinog::AbstractMatrix;
            rows::Optional{Integer} = nothing,
            cols::Optional{Integer} = nothing,
            α::Real = 360,
            α₀::Real = 0,
            ν::Real = 1,
            kwargs...
        )
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
            $f(sinog, xs, ys, ϕ₀..ϕ₁; kwargs...)
        end
    end
end


include("fbp_default.jl")
include("fbp_square.jl")
include("fbpa.jl")
include("fbpa_square.jl")


function iradon(
    sinog::AbstractMatrix,
    alg::AbstractIRadonAlgorithm = FBP();
    kwargs...
)
    alg(sinog; kwargs...)
end

function iradon(
    sinog::AbstractMatrix,
    xs::AbstractVector,
    ys::AbstractVector,
    alg::AbstractIRadonAlgorithm = FBP();
    kwargs...
)
    alg(sinog, xs, ys; kwargs...)
end

function iradon(
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
