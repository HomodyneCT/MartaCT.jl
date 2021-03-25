module CTImageData

export AbstractCTData, CTData

using ...Monads
using ...CTImages
import ...CTImages: ctimage, ctsinogram, cttomogram
import ...Utils: _atype
import Base: similar, eltype


abstract type AbstractCTData end


# Container type of CTData
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
mutable struct CTData{T<:Number,M<:AbstractMatrix{T}} <: AbstractCTData
    image::Maybe{CTImage{T,M}}
    sinog::Maybe{CTSinogram{T,M}}
    tomog::Maybe{CTTomogram{T,M}}
end


_atype(::Type{CTData{T,M}}) where {T,M} = M


function CTData{T,M}(gsd::AbstractCTData) where {T,M<:AbstractArray{T}}
    CTData{T,M}(
        mmap(CTImage{T,M}, gsd.image),
        mmap(CTSinogram{T,M}, gsd.sinog),
        mmap(CTTomogram{T,M}, gsd.tomog),
    )
end


CTData{T,M}() where {T,M<:AbstractArray{T}} =
    CTData{T,M}(nothing, nothing, nothing)


"""
    CTData(::Type{M}) where {M<:AbstractArray{<:Number}}

Construct an empty `CTData` object.
"""
CTData(::Type{M}) where {M<:AbstractArray{<:Number}} = CTData{eltype(M),M}()


"""
    CTData([::Type{T} = Float32]) where {T<:Number}

Construct an empty `CTData` object with eltype `T`.
"""
CTData(::Type{T} = Float32) where {T<:Number} = CTData(T,Matrix{T})


"""
    CTData(::AbstractCTImage)

Construct a `CTData` object from its constituents.
"""
CTData(img::CTImage) =
    CTData{eltype(img),_atype(img)}(Some(img), nothing, nothing)
CTData(img::CTSinogram) =
    CTData{eltype(img),_atype(img)}(nothing, Some(img), nothing)
CTData(img::CTTomogram) =
    CTData{eltype(img),_atype(img)}(nothing, nothing, Some(img))


function CTData{T,M}(gsd::AbstractCTData, img::CTImage) where {T,M}
    CTData{T,M}(
        Some(img),
        gsd.sinog,
        gsd.tomog,
    )
end


CTData{T,M}(gsd::AbstractCTData, img::CTSinogram) where {T,M} =
    CTData{T,M}(
        gsd.image,
        Some(img),
        gsd.tomog,
    )


CTData{T,M}(gsd::AbstractCTData, img::CTTomogram) where {T,M} =
    CTData{T,M}(
        gsd.image,
        gsd.sinog,
        Some(img),
    )


function CTData{T,M}(
    gsd::AbstractCTData,
    img::AbstractCTImage,
    rest::AbstractCTImage...
) where {T,M}
    CTData{T,M}(CTData{T,M}(gsd, img), rest...)
end

CTData(gsd::AbstractCTData, imgs::AbstractCTImage...) =
    CTData{eltype(gsd),_atype(gsd)}(gsd, imgs...)

function CTData(img::AbstractCTImage, imgs::AbstractCTImage...)
    CTData{eltype(img),_atype(img)}(CTData(img), imgs...)
end


similar(gsd::AbstractCTData, img::AbstractCTImage, rest::AbstractCTImage...) =
    typeof(gsd)(gsd, img, rest...)
similar(gsd::AbstractCTData) = typeof(gsd)()

end # module
