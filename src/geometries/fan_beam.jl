abstract type AbstractFanBeamGeometry <: AbstractGeometry end

struct FanBeamGeometry{T <: Real,CT <: AbstractTomograph} <:
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

@inline function getproperty(g::FanBeamGeometry, s::Symbol)
    s ≡ :width && return g.cols
    s ≡ :height && return g.rows
    s ≡ :nphi && return getfield(g, :nϕ)
    s ≡ :alpha && return getfield(g, :α)
    s ≡ :alpha0 && return getfield(g, :α₀)
    s ≡ :D1 && return getfield(g, :D′)
    s ≡ :gamma && return getfield(g, :γ)
    s ≡ :dx && return getfield(g, :δ)
    getfield(g, s)
end

eltype(::Type{G}) where {G <: FanBeamGeometry{T}} where {T} = T

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
    FanBeamGeometry([T=Float32]; <keyword arguments>) where {T}

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
    nϕ::Optional{Integer} = nothing,
    nphi::Optional{Integer} = nothing,
    D::Real = 500,
    D′::Optional{Real} = nothing,
    D1::Optional{Real} = nothing,
    γ::Optional{Real} = nothing,
    gamma::Optional{Real} = nothing,
    nd::Optional{Integer} = nothing,
    rows::Optional{Integer} = nothing,
    cols::Optional{Integer} = nothing,
    width::Optional{Integer} = nothing,
    height::Optional{Integer} = nothing,
    δ::Real = one(T),
    dx::Real = one(T),
    α::Real = 360,
    alpha::Real = 360,
    α₀::Real = zero(T),
    alpha0::Real = zero(T),
    center::Optional{Real} = nothing,
) where {T<:Real}
    nϕ = maybe(nϕ, nphi)
    D′ = maybe(D′, D1)
    γ = maybe(γ, gamma)
    δ = maybe(δ, dx)
    α = maybe(α, alpha)
    α₀ = maybe(α₀, alpha0)
    height = maybe(width, height)
    width = maybe(height, width)
    rows = maybe(512, maybe(height, rows))
    cols = maybe(rows, maybe(width, cols))
    nd = isnothing(nd) ? round(Int, hypot(rows, cols)) : nd
    nϕ = isnothing(nϕ) ?
        2 * (round(Int, (rows + 1) * (cols + 1) / nd) ÷ 2) + 1 : nϕ
    γ′::T = maybe(deg2rad, maybe(x->nd * δ / x, one(T), D′), γ)
    D′′::T = maybe(nd * δ / γ′, D′) # changed γ to γ′ here, should be tested
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


function show(io::IO, g::AbstractFanBeamGeometry)
    print(
    io,
    """Fan beam geometry:
    - tomograph: $(tomograph(g))
    - number of angles: $(num_proj(g))
    - number of detectors: $(num_det(g))
    - tomogram size (W×H): $(num_cols(g)) × $(num_rows(g))
    - focal spot to ISO: $(f2iso(g)) mm
    - fan beam angle: $(rad2deg(fan_angle(g)))°
    - cell size: $(cell_size(g)) mm
    - scan angle: $(scan_angle(g))°
    - scan starting angle: $(start_angle(g))°
    - channel spacing: $(maybe(x -> "$x°", "N/A", channel_spacing(g)))
    - center channel: $(center_channel(g))""",
    )
end
