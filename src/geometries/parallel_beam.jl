abstract type AbstractParallelBeamGeometry <: AbstractGeometry end

struct ParallelBeamGeometry{T <: Real,CT <: AbstractTomograph} <: AbstractParallelBeamGeometry
    ct::CT
    nϕ::Int
    nd::Int
    rows::Int
    cols::Int
    α::T
    α₀::T
    center::T
end

@inline function getproperty(g::ParallelBeamGeometry, s::Symbol)
    s ≡ :width && return g.cols
    s ≡ :height && return g.rows
    s ≡ :nphi && return getfield(g, :nϕ)
    s ≡ :alpha && return getfield(g, :α)
    s ≡ :alpha0 && return getfield(g, :α₀)
    getfield(g, s)
end

eltype(::Type{G}) where {G <: ParallelBeamGeometry{T}} where {T} = T

tomograph(x::ParallelBeamGeometry) = x.ct
num_proj(x::ParallelBeamGeometry) = x.nϕ
num_det(x::ParallelBeamGeometry) = x.nd
num_rows(x::ParallelBeamGeometry) = x.rows
num_cols(x::ParallelBeamGeometry) = x.cols
scan_angle(x::ParallelBeamGeometry) = x.α
start_angle(x::ParallelBeamGeometry) = x.α₀
center_channel(x::ParallelBeamGeometry) = x.center


"""
    ParallelBeamGeometry([T=Float32]; <keyword arguments>) where {T}

Construct geometry for the simulation.

# Arguments
- `nϕ`: number of projections.
- `nd=nothing`: number of detectors. If not specified, same as `nϕ`.
- `rows=nothing`: number of rows in the reconstructed image. If not specified, same as `nd`.
- `cols=nothing`: number of columns in the reconstructed image. If not
  specified, same as `rows`.
- `width=nothing`: alias for `cols`.
- `height=nothing`: alias for `rows`.
- `α=360`: scan angle in degrees.
- `α₀=0`: starting scan angle in degrees.
- `center=nothing`: virtual center channel, defaults to `(nd-1)/2`.
"""
function ParallelBeamGeometry(
    ::Type{T} = Float32,
    ct::AbstractTomograph = DefaultTomograph();
    nϕ::Optional{Integer} = nothing,
    nphi::Optional{Integer} = nothing,
    nd::Optional{Integer} = nothing,
    rows::Optional{Integer} = nothing,
    cols::Optional{Integer} = nothing,
    width::Optional{Integer} = nothing,
    height::Optional{Integer} = nothing,
    α::Real = 360,
    alpha::Optional{Real} = nothing,
    α₀::Real = 0,
    alpha0::Optional{Real} = nothing,
    center::Optional{Real} = nothing,
) where {T<:Real}
    nϕ = maybe(nϕ, nphi)
    α = maybe(α, alpha)
    α₀ = maybe(α₀, alpha0)
    height = maybe(width, height)
    width = maybe(height, width)
    rows = maybe(512, maybe(height, rows))
    cols = maybe(rows, maybe(width, cols))
    nd = isnothing(nd) ? round(Int, hypot(rows, cols)) : nd
    nϕ = isnothing(nϕ) ?
        2 * (round(Int, (rows + 1) * (cols + 1) / nd) ÷ 2) + 1 : nϕ
    center′::T = maybe((nd - 1) / 2, center)
    ParallelBeamGeometry{T,typeof(ct)}(
        ct,
        nϕ,
        nd,
        rows,
        cols,
        α,
        α₀,
        center′,
    )
end


ParallelBeamGeometry(ct::CT; kwargs...) where {CT <: AbstractTomograph} =
    ParallelBeamGeometry(Float32, ct; kwargs...)


show(io::IO, g::AbstractParallelBeamGeometry) = print(
    io,
    """Parallel beam geometry:
- tomograph: $(g.ct)
- number of projections: $(num_proj(g))
- number of detectors: $(num_det(g))
- tomogram size (W×H): $(num_cols(g)) × $(num_rows(g))
- scan angle: $(scan_angle(g))°
- scan starting angle: $(start_angle(g))°
- center channel: $(center_channel(g))""",
)
