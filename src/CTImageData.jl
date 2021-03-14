module CTImageData

export AbstractCTData, CTData

using ..Monads
using ..CTImages
import ..CTImages: ctimage, ctsinogram, cttomogram
import ..Marta: _atype
import Base: similar, eltype


abstract type AbstractCTData{M<:AbstractArray{<:Number}} end


# Container type of CTData
_atype(::Type{T}) where {T<:AbstractCTData{M}} where {M} = M
_atype(x::AbstractCTData) = _atype(typeof(x))

eltype(::Type{T}) where {T<:AbstractCTData} = eltype(_atype(T))

ctimage(gsd::AbstractCTData) = gsd.image
ctsinogram(gsd::AbstractCTData) = gsd.sinog
cttomogram(gsd::AbstractCTData) = gsd.tomog


"""
    mutable struct CTData{M}

Hold data for reconstruction.
`M` should be a subtype of `AbstractArray`, for example
`Matrix{Float32}` or `AFMatrix{Float64}`.
By default, the library will use CPU objects rather than GPU ones.
"""
mutable struct CTData{M} <: AbstractCTData{M}
    image::Maybe{CTImage{M}}
    sinog::Maybe{CTSinogram{M}}
    tomog::Maybe{CTTomogram{M}}
end


function CTData{M}(gsd::AbstractCTData) where {M<:AbstractArray{<:Number}}
    CTData{M}(
        mmap(CTImage{M}, gsd.image),
        mmap(CTSinogram{M}, gsd.sinog),
        mmap(CTTomogram{M}, gsd.tomog),
    )
end


CTData{M}() where {M<:AbstractArray{<:Number}} =
    CTData{M}(nothing, nothing, nothing)


"""
    CTData(::Type{M}) where {M<:AbstractArray{<:Number}}

Construct an empty `CTData` object.
"""
CTData(::Type{M}) where {M<:AbstractArray{<:Number}} =
    CTData{M}()


"""
    CTData([::Type{T} = Float32]) where {T<:Number}

Construct an empty `CTData` object with eltype `T`.
"""
CTData(::Type{T} = Float32) where {T<:Number} =
    CTData(Matrix{T})


"""
    CTData(::AbstractCTImage{M}) where {M<:AbstractArray{<:Number}}

Construct a `CTData` object from its constituents.
"""
CTData(img::CTImage) = CTData{_atype(img)}(Some(img), nothing, nothing)
CTData(img::CTSinogram) = CTData{_atype(img)}(nothing, Some(img), nothing)
CTData(img::CTTomogram) = CTData{_atype(img)}(nothing, nothing, Some(img))


function CTData{M}(gsd::AbstractCTData, img::CTImage) where {M}
    CTData{M}(
        Some(img),
        gsd.sinog,
        gsd.tomog,
    )
end


CTData{M}(gsd::AbstractCTData, img::CTSinogram) where {M} =
    CTData{M}(
        gsd.image,
        Some(img),
        gsd.tomog,
    )


CTData{M}(gsd::AbstractCTData, img::CTTomogram) where {M} =
    CTData{M}(
        gsd.image,
        gsd.sinog,
        Some(img),
    )


CTData{M}(gsd::AbstractCTData, img::AbstractCTImage, rest::AbstractCTImage...) where {M} =
    CTData{M}(CTData{M}(gsd, img), rest...)

CTData(gsd::AbstractCTData, imgs::AbstractCTImage...) =
    CTData{_atype(gsd)}(gsd, imgs...)

function CTData(img::AbstractCTImage, imgs::AbstractCTImage...)
    CTData{_atype(img)}(CTData(img), imgs...)
end


similar(gsd::AbstractCTData, imgs::AbstractCTImage...) = typeof(gsd)(gsd, imgs...)
