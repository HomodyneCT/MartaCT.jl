abstract type AbstractConeBeamGeometry <: AbstractGeometry end

struct ConeBeamFlatGeometry{T <: Real,Q <: AbstractQuantity{T},CT <: AbstractTomograph} <: AbstractConeBeamGeometry
    ct::CT
    isize::NTuple{3,Int}
    osize::NTuple{3,Int}
    DSO::Q
    δ::Q
    center::T
end


eltype(::Type{G}) where {G <: ConeBeamFlatGeometry{T}} where {T} = T

tomograph(x::ConeBeamFlatGeometry) = x.ct
num_proj(x::ConeBeamFlatGeometry) = x.isize[3]
num_det(x::ConeBeamFlatGeometry) = x.isize[begin:end-1]
f2iso(x::ConeBeamFlatGeometry) = x.DSO
cell_size(x::ConeBeamFlatGeometry) = x.δ
channel_spacing(x::ConeBeamFlatGeometry) = tomograph(x).channel_spacing
center_channel(x::ConeBeamFlatGeometry) = x.center


function ConeBeamFlatGeometry(
    ::Type{T} = Float32,
    ct::AbstractTomograph = DefaultTomograph();
    isize::NTuple{3,Int},
    osize::NTuple{3,Int} = (512,512,512),
    DSO::Optional{AbstractQuantity} = nothing,
    R::Optional{AbstractQuantity} = nothing,
    δ::Optional{AbstractQuantity} = nothing,
    dx::Optional{AbstractQuantity} = nothing,
    center::Optional{Real} = nothing,
) where {T<:Real}
    DSO = maybe(R, DSO)
    δ = uconvert(unit(DSO), maybe(dx, δ))
    center′::T = maybe((nd - 1) / 2, center)
    ConeBeamFlatGeometry{T,eltype(DSO),typeof(ct)}(
        ct,
        isize,
        osize,
        DSO,
        δ,
        center′
    )
end


ConeBeamFlatGeometry(ct::AbstractTomograph; kwargs...) =
    ConeBeamFlatGeometry(Float32, ct; kwargs...)
