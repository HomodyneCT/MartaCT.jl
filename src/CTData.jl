module CTData

export AbstractCTData, GrayScaleData


using ..Monads
using ..CTImages
import ..CTImages: ctimage, ctsinogram, cttomogram
import ..Marta: _atype
import Base: similar


abstract type AbstractCTData{M<:AbstractArray{<:Number}} end


# Container type of CTData
_atype(::Type{<:AbstractCTData{M}}) where {M} = M
_atype(x::AbstractCTData) = _atype(typeof(x))


ctimage(gsd::AbstractCTData) = gsd.image
ctsinogram(gsd::AbstractCTData) = gsd.sinog
cttomogram(gsd::AbstractCTData) = gsd.tomog


"""
    struct GrayScaleData{M}

Hold data for reconstruction.
`M` should be a subtype of `AbstractArray`, for example
`Matrix{Float32}` or `AFMatrix{Float64}`.
By default, the library will use CPU objects rather than GPU ones.
"""
struct GrayScaleData{M} <: AbstractCTData{M}
    image::Maybe{CTImage{M}}
    sinog::Maybe{CTSinogram{M}}
    tomog::Maybe{CTTomogram{M}}
end


function GrayScaleData{M}(gsd::AbstractCTData) where {M<:AbstractArray{<:Number}}
    GrayScaleData{M}(
        mmap(CTImage{M}, gsd.image),
        mmap(CTSinogram{M}, gsd.sinog),
        mmap(CTTomogram{M}, gsd.tomog),
    )
end


GrayScaleData{M}() where {M<:AbstractArray{<:Number}} =
    GrayScaleData{M}(nothing, nothing, nothing)


"""
    GrayScaleData(::Type{M}) where {M<:AbstractArray{<:Number}}

Construct an empty `GrayScaleData` object.
"""
GrayScaleData(::Type{M}) where {M<:AbstractArray{<:Number}} =
    GrayScaleData{M}()


"""
    GrayScaleData([::Type{T} = Float32]) where {T<:Number}

Construct an empty `GrayScaleData` object with eltype `T`.
"""
GrayScaleData(::Type{T} = Float32) where {T<:Number} =
    GrayScaleData(Matrix{T})


"""
    GrayScaleData(::AbstractCTImage{M}) where {M<:AbstractArray{<:Number}}

Construct a `GrayScaleData` object from its constituents.
"""
GrayScaleData(img::CTImage) = GrayScaleData{_atype(img)}(Some(img), nothing, nothing)
GrayScaleData(img::CTSinogram) = GrayScaleData{_atype(img)}(nothing, Some(img), nothing)
GrayScaleData(img::CTTomogram) = GrayScaleData{_atype(img)}(nothing, nothing, Some(img))


function GrayScaleData{M}(gsd::AbstractCTData, img::CTImage) where {M}
    GrayScaleData{M}(
        Some(img),
        gsd.sinog,
        gsd.tomog,
    )
end


GrayScaleData{M}(gsd::AbstractCTData, img::CTSinogram) where {M} =
    GrayScaleData{M}(
        gsd.image,
        Some(img),
        gsd.tomog,
    )


GrayScaleData{M}(gsd::AbstractCTData, img::CTTomogram) where {M} =
    GrayScaleData{M}(
        gsd.image,
        gsd.sinog,
        Some(img),
    )


GrayScaleData{M}(gsd::AbstractCTData, img::AbstractCTImage, rest::AbstractCTImage...) where {M} =
    GrayScaleData{M}(GrayScaleData{M}(gsd, img), rest...)

GrayScaleData(gsd::AbstractCTData, imgs::AbstractCTImage...) =
    GrayScaleData{_atype(gsd)}(gsd, imgs...)

function GrayScaleData(img::AbstractCTImage, imgs::AbstractCTImage...)
    GrayScaleData{_atype(img)}(GrayScaleData(img), imgs...)
end


similar(gsd::AbstractCTData, imgs::AbstractCTImage...) = typeof(gsd)(gsd, imgs...)

end # module
