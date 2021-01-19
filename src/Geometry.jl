module Geometry

Base.Experimental.@optlevel 1

export AbstractGeometry,
    AbstractParallelBeamGeometry, AbstractFanBeamGeometry
export ParallelBeamGeometry, FanBeamGeometry
export AbstractTomograph, DefaultTomograph
export num_proj, num_det, f2iso, f2det, fan_angle, cell_size
export scan_angle, start_angle, center_channel
export num_rows, num_cols, tomograph, channel_spacing

import Base: show, getproperty

using ..Monads
import ..CTIO: yaml_repr, struct2dict
import ..Marta: datatype


abstract type AbstractGeometry  end

datatype(x::AbstractGeometry) = datatype(typeof(x))

function getproperty(g::AbstractGeometry, s::Symbol)
    s ≡ :width && return g.cols
    s ≡ :height && return g.rows
    getfield(g, s)
end


abstract type AbstractTomograph end

struct DefaultTomograph <: AbstractTomograph end


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

datatype(::Type{<:ParallelBeamGeometry{T}}) where {T} = T

yaml_repr(g::ParallelBeamGeometry) = struct2dict(g)

tomograph(x::ParallelBeamGeometry) = x.ct
num_proj(x::ParallelBeamGeometry) = x.nϕ
num_det(x::ParallelBeamGeometry) = x.nd
num_rows(x::ParallelBeamGeometry) = x.rows
num_cols(x::ParallelBeamGeometry) = x.cols
scan_angle(x::ParallelBeamGeometry) = x.α
start_angle(x::ParallelBeamGeometry) = x.α₀
center_channel(x::ParallelBeamGeometry) = x.center


"""
    ParallelBeamGeometry([T=Float32]; <keyword arguments>) where {T<:Real}

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
    nϕ::Int = 1024,
    nd::Optional{Int} = nothing,
    rows::Optional{Int} = nothing,
    cols::Optional{Int} = nothing,
    width::Optional{Int} = nothing,
    height::Optional{Int} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    center::Optional{<:Real} = nothing,
) where {T <: Real}
    nd = maybe(nϕ, nd)
    rows = maybe(nd, maybe(rows, height))
    cols = maybe(rows, maybe(cols, width))
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



abstract type AbstractFanBeamGeometry <: AbstractGeometry end

struct FanBeamGeometry{T,CT <: AbstractTomograph} <:
       AbstractFanBeamGeometry
    ct::CT
    nϕ::Int
    nd::Int
    rows::Int
    cols::Int
    D::T
    D′::T
    γ::T
    δ::T
    α::T
    α₀::T
    center::T
end


datatype(::Type{<:FanBeamGeometry{T}}) where {T} = T

yaml_repr(g::FanBeamGeometry) = struct2dict(g)

tomograph(x::FanBeamGeometry) = x.ct
num_proj(x::FanBeamGeometry) = x.nϕ
num_det(x::FanBeamGeometry) = x.nd
num_rows(x::FanBeamGeometry) = x.rows
num_cols(x::FanBeamGeometry) = x.cols
f2iso(x::FanBeamGeometry) = x.D
f2det(x::FanBeamGeometry) = x.D′
fan_angle(x::FanBeamGeometry) = x.γ
cell_size(x::FanBeamGeometry) = x.δ
scan_angle(x::FanBeamGeometry) = x.α
start_angle(x::FanBeamGeometry) = x.α₀
channel_spacing(x::FanBeamGeometry) = tomograph(x).channel_spacing
center_channel(x::FanBeamGeometry) = x.center


"""
    FanBeamGeometry([T=Float32]; <keyword arguments>) where {T<:Real}

Create a new FanBeamGeometry object with given parameters.

It is assumed that the detectors are arranged on a circle arc,
not flat.

# Arguments
- `nϕ: number of scan angles.
- `nd=nothing`: number of detectors.
- `rows=nothing`: number of rows in the reconstructed image.
- `cols=nothing`: number of columns in the reconstructed image.
- `width=nothing`: alias for `cols`.
- `height=nothing`: alias for `rows`.
- `D`: focal spot to ISO distance.
- `D′=nothing`: focal spot to detectors distance; if not specified
    and also `γ` is not specified, is defined so that the total
    fan angle is 1 radian.
- `γ=nothing`: total fan angle in degrees.
- `δ=1`: detectors spacing (cell size).
- `α=360`: scan angle in degrees.
- `α₀=0`: scan starting angle in degrees.
- `center=nothing`: virtual center channel, defaults to `(nd-1)/2`.
"""
function FanBeamGeometry(
    ::Type{T} = Float32,
    ct::AbstractTomograph = DefaultTomograph();
    nϕ::Int = 1024,
    D::Real,
    D′::Optional{<:Real} = nothing,
    γ::Optional{<:Real} = nothing,
    nd::Optional{Int} = nothing,
    rows::Optional{Int} = nothing,
    cols::Optional{Int} = nothing,
    width::Optional{Int} = nothing,
    height::Optional{Int} = nothing,
    δ::Real = 1,
    α::Real = 360,
    α₀::Real = 0,
    center::Optional{<:Real} = nothing,
) where {T <: Real}
    nd = maybe(nϕ, nd)
    rows = maybe(nd, maybe(rows, height))
    cols = maybe(rows, maybe(cols, width))
    γ′::T = maybe(
        deg2rad,
        maybe(1, D′) do x
            nd * δ / x
        end,
        γ
    )
    D′′::T = maybe(nd * δ / γ, D′)
    center′::T = maybe((nd - 1) / 2, center)
    FanBeamGeometry{T,typeof(ct)}(
        ct,
        nϕ,
        nd,
        rows,
        cols,
        D,
        D′′,
        γ′,
        δ,
        α,
        α₀,
        center′
    )
end


FanBeamGeometry(ct::AbstractTomograph; kwargs...) =
    FanBeamGeometry(Float32, ct; kwargs...)


show(io::IO, g::AbstractFanBeamGeometry) = print(
    io,
    """Fan beam geometry:
- tomograph:
    $(tomograph(g))
- number of angles: $(num_proj(g))
- number of detectors: $(num_det(g))
- tomogram size (W×H): $(num_cols(g)) × $(num_rows(g))
- focal spot to ISO: $(f2iso(g)) mm
- fan beam angle: $(rad2deg(fan_angle(g)))°
- cell size: $(cell_size(g)) mm
- scan angle: $(scan_angle(g))°
- scan starting angle: $(start_angle(g))°
- channel spacing: $(channel_spacing(g))°
- center channel: $(center_channel(g))""",
)



function ParallelBeamGeometry(fbg::AbstractFanBeamGeometry)
    ParallelBeamGeometry(
        datatype(fbg),
        fbg.ct;
        fbg.nϕ,
        fbg.nd,
        fbg.rows,
        fbg.cols,
        fbg.α,
        fbg.α₀,
        fbg.center,
    )
end


function FanBeamGeometry(
    pbg::AbstractParallelBeamGeometry,
    ct::AbstractTomograph = DefaultTomograph();
    kwargs...,
)
    FanBeamGeometry(
        datatype(pbg),
        ct;
        pbg.nϕ,
        pbg.nd,
        pbg.rows,
        pbg.cols,
        pbg.α,
        pbg.α₀,
        pbg.center,
        kwargs...,
    )
end

end # module
