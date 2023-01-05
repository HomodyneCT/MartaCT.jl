module CTImages

export rescale!, rescale, rotate
export polar2cart
# export AbstractCTImage, CTImage, CTSinogram, CTTomogram
# export CTImageOrTomog, CTImageMat, CTSinogMat, CTTomogMat
# export ctimage, ctsinogram, cttomogram

using IntervalSets, SimpleTraits
import ..Monads: mjoin
using MartaCT.Utils: linspace, half, _compute_radius, _compute_angle
import MartaCT.Utils: _atype
using MartaCT.Applicative, ..Monads, ..Interpolation
using MartaCT.Geometry:
    AbstractParallelBeamGeometry, AbstractFanBeamGeometry, ParallelBeamGeometry
import Base: size, convert, similar, eltype, ndims, getindex, setindex!, IndexStyle
using Base: @propagate_inbounds


function ctimage end
function ctsinogram end
function cttomogram end


abstract type AbstractCTImage{T,N} <: AbstractArray{T,N} end

@traitimpl Monad{AbstractCTImage}
mjoin(img::AbstractCTImage) = img.data

# struct CTImageStyle <: Broadcast.AbstractArrayStyle{2} end
# Base.BroadcastStyle(::Type{<:AbstractCTImage}) = CTImageStyle()
# CTImageStyle(::Val{2}) = CTImageStyle()
# CTImageStyle(::Val{N}) where N = Broadcast.DefaultArrayStyle{N}()

_atype(img::M) where {M <: AbstractCTImage} = _atype(M)
eltype(::Type{M}) where {M<:AbstractCTImage{T}} where T = T
ndims(::Type{M}) where {M<:AbstractCTImage{T,N}} where {T,N} = ndims(_atype(M))
size(img::AbstractCTImage) = mbind(size, img)
size(img::AbstractCTImage, dim::Int) = size(mjoin(img), dim)
@propagate_inbounds function getindex(img::AbstractCTImage, i::Int)
    @boundscheck checkbounds(img, i)
    @inbounds getindex(mjoin(img), i)
end
@propagate_inbounds function getindex(
    img::AbstractCTImage, I::Vararg{Int,N}
) where N
    @boundscheck checkbounds(img, I...)
    @inbounds getindex(mjoin(img), I...)
end
@propagate_inbounds function setindex!(img::AbstractCTImage, v, i::Int)
    @boundscheck checkbounds(img, i)
    @inbounds setindex!(mjoin(img), v, i)
end
@propagate_inbounds function setindex!(
    img::AbstractCTImage, v, I::Vararg{Int,N}
) where N
    @boundscheck checkbounds(img, I...)
    @inbounds setindex!(mjoin(img), v, I...)
end


const ctnames = (image = :CTImage, sinog = :CTSinogram, tomog = :CTTomogram)
const ctfn = (image = :ctimage, sinog = :ctsinogram, tomog = :cttomogram)


for nm in ctnames
    @eval begin
        struct $nm{T,N,M<:AbstractArray{T,N}} <: AbstractCTImage{T,N}
            data::M
        end
        @inline _atype(::Type{$nm{T,N,M}}) where {T,N,M} = M
        IndexStyle(::Type{T}) where T<:$nm = IndexStyle(_atype(T))
        @inline function $nm{T,N,M}(
            ::UndefInitializer,
            dims::Vararg{Union{Integer,AbstractUnitRange},N}
        ) where {T,N,M}
            $nm(M(undef, dims...))
        end
        @inline $nm{T,N,M}(img::$nm) where {T,N,M} = $nm{T,N,M}(mjoin(img))
        @inline $nm{T}(img::M) where {T,M<:AbstractArray{T,N}} where N = $nm{T,N,M}(img)
        @inline $nm(img::$nm) = img
        @inline convert(::Type{M}, img::AbstractArray) where {M<:$nm} = mreturn(M, img)
        @inline convert(::Type{M}, img::$nm) where {M<:$nm} = mreturn($nm, convert(_atype(M), mjoin(img)))
        @inline similar(img::M) where {M<:$nm} = mreturn(M, similar(mjoin(img)))
        @inline function similar(
            img::M,
            ::Type{S},
            dims::Tuple{Vararg{Int,N}}
        ) where {S,N,M<:$nm}
            mreturn(M, similar(mjoin(img), S, dims))
        end
    end
end

# CTImage(::Union{CTSinogram,CTTomogram}) =
#     @error "Cannot construct a CT image from a sinogram or a tomogram"
# CTSinogram(::Union{CTImage,CTTomogram}) =
#     @error "Cannot construct a CT sinogram from an image or a tomogram"
# CTTomogram(::Union{CTImage,CTSinogram}) =
#     @error "Cannot construct a CT tomogram from an image or a sinogram"

const CTImageOrTomog{T,N,M} = Union{CTImage{T,N,M},CTTomogram{T,N,M}}
const CTImageMat{T} = CTImage{T,2,Matrix{T}}
const CTSinogMat{T} = CTSinogram{T,2,Matrix{T}}
const CTTomogMat{T} = CTTomogram{T,2,Matrix{T}}

"""
    rotate(mat::AbstractMatrix, α::Real; <keyword arguments>)

Rotate matrix `mat` about the center of angle `α` given in degrees.
If `rows` and `cols` are not given, the rotated matrix has the same
dimensions of the original matrix.

# Arguments
- `mat`: matrix to rotate.
- `α`: angle in degrees.
- `rows=nothing`: number of rows of the rotated matrix.
- `cols=nothing`: number of columns of the rotated matrix.
- `interpolation`: interpolation strategy. By default is
    `BilinearInterpolation`.
"""
function rotate(
    mat::AbstractMatrix{T},
    α::Real;
    rows::Optional{Integer} = nothing,
    cols::Optional{Integer} = nothing,
    background::Optional{Real} = nothing,
    interpolation::Optional{Interp} = nothing,
) where {T,Interp <: AbstractInterp2DOrNone}
    Base.depwarn("`rotate` will be removed in a future release, please use the `Images` package functionality", :rotate)
    orows, ocols = size(mat)
    rows = maybe(orows, rows)
    cols = maybe(ocols, cols)
    sϕ, cϕ = sincos(deg2rad(α))
    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(mat)
    x₀::T = (cols + 1) / 2
    y₀::T = (rows + 1) / 2
    x′₀::T = (ocols + 1) / 2
    y′₀::T = (orows + 1) / 2
    z::T = maybe(zero(T), background)
    rmat = similar(mat, rows, cols)
    fill!(rmat, z)
    for ix ∈ 1:rows, iy ∈ 1:cols
        x = T(ix) - x₀
        y = T(iy) - y₀
        x′ = x * cϕ + y * sϕ + x′₀
        y′ = y * cϕ - x * sϕ + y′₀
        if x′ ∈ 1..ocols && y′ ∈ 1..orows
            rmat[iy, ix] = convert(T, interp[y′, x′])
        end
    end
    rmat
end


"""
    rescale!(x::AbstractArray;
        interval=nothing, calibration=nothing, window=nothing)

Rescale array x to the interval specified by `interval`.

If `calibration` is not `nothing`, then rescaling is done with
reference to the values given by `calibration`. In other words,
minimum and maximum are assumed to be the values specified by
`calibration`.

See also: [`rescale`](@ref)
"""
function rescale!(
    img::AbstractArray{T};
    interval::ClosedInterval = zero(T)..one(T),
    calibration::Optional{ClosedInterval} = nothing,
    window::Optional{ClosedInterval} = nothing,
) where {T}
    Base.depwarn("`rescale` will be removed in a future release, please use the `Images` package functionality", :rescale)
    if isnothing(calibration)
        calibration = ClosedInterval(extrema(img)...)
    end
    a, b = T.(endpoints(interval))
    m, M = T.(endpoints(calibration))
    if m == M
        @warn "Cannot calibrate image as calibration values ($calibration) are equal, image unchanged"
        return img
    end
    f = (b - a) / (M - m)
    if m != zero(T)
        img .-= m
    end
    img .*= f
    if a != zero(T)
        img .+= a
    end
    isnothing(window) && return img
    a, b = T.(endpoints(window))
    map!(img, img) do x
        x < a && return a
        x > b && return b
        return x
    end
    img
end



"""
    rescale(x::AbstractArray, slope::Number, intercept::Number; window)

Linear rescaling of `x` as `x * slope + intercept`.

See also: [`rescale!`](@ref)
"""
function rescale(
    image::AbstractArray{T},
    slope::Number,
    intercept::Number;
    window::Optional{ClosedInterval} = nothing,
) where {T}
    Base.depwarn("`rescale` will be removed in a future release, please use the `Images` package functionality", :rescale)
    res = @. image * slope + intercept
    isnothing(window) && return res
    a, b = T.(endpoints(window))
    map!(res, res) do x
        x < a && return a
        x > b && return b
        return x
    end
end



"""
    rescale!(x::AbstractArray, slope::Number, intercept::Number)

In place linear rescaling of `x`.

See also: [`rescale`](@ref)
"""
function rescale!(
    img::AbstractArray{T},
    slope::Number,
    intercept::Number;
    window::Optional{ClosedInterval} = nothing,
) where {T}
    Base.depwarn("`rescale` will be removed in a future release, please use the `Images` package functionality", :rescale)
    @. img = slope * img + intercept
    isnothing(window) && return img
    a, b = T.(endpoints(window))
    map!(img, img) do x
        x < a && return a
        x > b && return b
        return x
    end
    img
end



"""
    rescale(x::AbstractArray;
        interval=nothing, calibration=nothing, window=nothing)

Rescale array x to the interval specified by `interval`.

If `calibration` is not `nothing`, then rescaling is done with
reference to the values given by `calibration`. In other words,
minimum and maximum are assumed to be the values specified by
`calibration`.

See also: [`rescale!`](@ref)
"""
function rescale(
    x::AbstractArray{T};
    interval = zero(T)..one(T),
    calibration = nothing,
    window = nothing,
) where {T}
    rescale!(deepcopy(x); interval, calibration, window)
end


function polar2cart(
    mp::AbstractMatrix{T},
    xs::AbstractVector,
    ys::AbstractVector;
    ν::Real = 1,
    scale::Optional{Real} = nothing,
    diag::Bool = false,
    background::Optional{Real} = nothing,
    transposed::Bool = false,
    interpolation::Optional{Interp} = nothing,
) where {T,Interp <: AbstractInterp2DOrNone}
    if transposed
        nr, nθ = size(mp)
        mp′ = similar(mp, nθ, nr)
        permutedims!(mp′, mp, (2,1))
        return polar2cart(mp′, xs, ys; background, interpolation)
    end
    ν = maybe(ν, scale)
    @assert ν > 0
    nθ, nr = size(mp)
    x₀, y₀ = half(xs), half(ys)
    κ::T = min(x₀, y₀)
    if diag
        ν *= hypot(x₀, y₀) / κ
    end
    δri::T = (nr-1) / (κ * ν)
    δθi::T = nθ / 2π
    xs′ = xs * δri
    ys′ = ys * δri
    nt1 = nθ + 1
    nthalf = nθ + 1//2
    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(mp)
    rows, cols = length(ys), length(xs)
    z::T = maybe(zero(T), background)
    mc = similar(mp, rows, cols)
    fill!(mc, z)
    Threads.@threads for ix ∈ axes(mc, 2)
        @inbounds @fastmath @simd for iy ∈ axes(mc, 1)
            x::T = xs′[ix]
            y::T = ys′[iy]
            r::T = _compute_radius(x, y)
            θ::T = _compute_angle(x, y, δθi)
            if 1 <= r <= nr && 1 <= θ <= nt1
                if θ >= nthalf
                    θ = one(T)
                elseif θ > nθ
                    θ = nθ
                end
                mc[iy,ix] = interp[θ, r]
            end
        end
    end
    mc
end


@inline function polar2cart(
    mp::AbstractMatrix{T};
    rows::Optional{Integer} = nothing,
    cols::Optional{Integer} = nothing,
    transposed::Bool = false,
    kwargs...
) where {T}
    mp = transposed ? permutedims(mp) : mp
    rows = isnothing(rows) ? maybe(2 * size(mp, 2), cols) : rows
    cols = maybe(rows, cols)
    y₀, x₀ = sincos(atan(rows, cols))
    xs = linspace(T, -x₀..x₀, cols)
    ys = linspace(T, -y₀..y₀, rows)
    polar2cart(mp, xs, ys; kwargs...)
end


function polar2cart(
    image::AbstractMatrix,
    geometry::AbstractParallelBeamGeometry;
    kwargs...
)
    polar2cart(image; geometry.rows, geometry.cols, kwargs...)
end

function polar2cart(
    image::AbstractMatrix,
    geometry::AbstractFanBeamGeometry;
    kwargs...
)
    polar2cart(image, ParallelBeamGeometry(geometry); kwargs...)
end

polar2cart(; kwargs...) = x -> polar2cart(x; kwargs...)
polar2cart(g::AbstractParallelBeamGeometry; kwargs...) =
    x -> polar2cart(x, g; kwargs...)

end # module
