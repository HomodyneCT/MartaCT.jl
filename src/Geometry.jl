module Geometry

Base.Experimental.@optlevel 1

export
    AbstractGeometry,
    AbstractParallelBeamGeometry,
    AbstractFanBeamGeometry,
    AbstractConeBeamGeometry
export ParallelBeamGeometry, FanBeamGeometry, ConeBeamFlatGeometry
export AbstractTomograph, DefaultTomograph
export num_proj, num_det, f2iso, f2det, fan_angle, cell_size
export scan_angle, start_angle, center_channel
export num_rows, num_cols, tomograph, channel_spacing

using ..Monads
import Base: show, getproperty, eltype
using Unitful: AbstractQuantity


abstract type AbstractTomograph end

function getproperty(ct::AbstractTomograph, s::Symbol)
    s ∈ fieldnames(typeof(ct)) && return getfield(ct, s)
    nothing
end

show(io::IO, ::T) where T <: AbstractTomograph = print(io, "$(nameof(T))")

struct DefaultTomograph <: AbstractTomograph end


const _geometry_names = (:ParallelBeamGeometry, :FanBeamGeometry, :ConeBeamFlatGeometry)

abstract type AbstractGeometry end

include("geometries/parallel_beam.jl")
include("geometries/fan_beam.jl")
include("geometries/cone_beam.jl")


function ParallelBeamGeometry(fbg::AbstractFanBeamGeometry)
    ParallelBeamGeometry(
        eltype(fbg),
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
        eltype(pbg),
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
