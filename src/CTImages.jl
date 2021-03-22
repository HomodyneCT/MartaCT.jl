module CTImages

export rescale!, rescale, rotate
export polar2cart
export AbstractCTImage, CTImage, CTSinogram, CTTomogram
export CTImageOrTomog, CTImageMat, CTSinogMat, CTTomogMat
export ctimage, ctsinogram, cttomogram

using IntervalSets, SimpleTraits
import ..Monads: mjoin
using ..Utils: linspace
import ..Utils: _atype
using ..Applicative, ..Monads, ..Interpolation
using ..Geometry:
    AbstractParallelBeamGeometry, AbstractFanBeamGeometry, ParallelBeamGeometry
import Base: size, convert, similar, eltype, getindex, setindex!, IndexStyle
using Base: @propagate_inbounds


function ctimage end
function ctsinogram end
function cttomogram end


abstract type AbstractCTImage{T<:Number} <: AbstractMatrix{T} end

@traitimpl Monad{AbstractCTImage}

mjoin(img::AbstractCTImage) = img.data

_atype(img::M) where {M <: AbstractCTImage} = _atype(M)
eltype(::Type{M}) where {M<:AbstractCTImage{T}} where T = T
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
        struct $nm{T,M<:AbstractMatrix{T}} <: AbstractCTImage{T}
            data::M
        end
        @inline _atype(::Type{$nm{T,M}}) where {T,M} = M
        IndexStyle(::Type{T}) where T<:$nm = IndexStyle(_atype(T))
        @inline function $nm{T,M}(
            ::UndefInitializer,
            dims::Vararg{Union{Integer,AbstractUnitRange}}
        ) where {T,M}
            $nm(M(undef, dims...))
        end
        @inline $nm{T,M}(img::$nm) where {T,M} = $nm{T,M}(mjoin(img))
        @inline $nm(img::$nm) = img
        @inline convert(::Type{M}, img::$nm) where {M<:$nm} =
            mreturn(M, convert(_atype(M), mjoin(img)))
        @inline function similar(
            ::Type{M},
            dims::Vararg{Union{Integer,AbstractUnitRange},N},
            ) where {M<:$nm,N}
            mreturn(M, similar(_atype(M), dims...))
        end
        @inline function similar(
            ::M,
            dims::Vararg{Union{Integer,AbstractUnitRange},N},
        ) where {M<:$nm,N}
            similar(M, dims...)
        end
    end
end

CTImage(::Union{CTSinogram,CTTomogram}) =
    @error "Cannot construct a CT image from a sinogram or a tomogram"
CTSinogram(::Union{CTImage,CTTomogram}) =
    @error "Cannot construct a CT sinogram from an image or a tomogram"
CTTomogram(::Union{CTImage,CTSinogram}) =
    @error "Cannot construct a CT tomogram from an image or a sinogram"

const CTImageOrTomog{T,M} = Union{CTImage{T,M},CTTomogram{T,M}}
const CTImageMat{T} = CTImage{T,Matrix{T}}
const CTSinogMat{T} = CTSinogram{T,Matrix{T}}
const CTTomogMat{T} = CTTomogram{T,Matrix{T}}

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
) where {T <: Number,Interp <: AbstractInterp2DOrNone}
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
            rmat[iy, ix] = interp(y′, x′)
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
    calibration::Optional{ClosedInterval{U}} = nothing,
    window::Optional{ClosedInterval{W}} = nothing,
) where {T <: Number,U <: Number,W <: Number}
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
    window::Optional{ClosedInterval{U}} = nothing,
) where {T <: Number,U <: Number}
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
    window::Optional{ClosedInterval{U}} = nothing,
) where {T <: Number,U <: Number}
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
    x::AbstractArray;
    interval = 0..1,
    calibration = nothing,
    window = nothing,
)
    rescale!(deepcopy(x); interval, calibration, window)
end


function polar2cart(
    mp::AbstractMatrix{T},
    xs::AbstractVector,
    ys::AbstractVector;
    background::Optional{Real} = nothing,
    transposed::Bool = false,
    interpolation::Optional{Interp} = nothing,
) where {T <: Real,Interp <: AbstractInterp2DOrNone}
    mp = transposed ? permutedims(mp) : mp
    nθ, nr = size(mp)
    Δr = nr - 1
    Δθ = (nθ - 1) / 2π
    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(mp)
    rows, cols = length(ys), length(xs)
    indices = Vector{NTuple{2,Int}}(undef, rows * cols)
    @inbounds for k ∈ eachindex(indices)
        indices[k] = (k - 1) ÷ rows + 1, (k - 1) % rows + 1
    end
    z::T = maybe(zero(T), background)
    mc = similar(mp, rows, cols)
    fill!(mc, z)
    @inline function _compute_radius(x, y)
        x == 0 && return abs(y)
        y == 0 && return abs(x)
        return √(x^2 + y^2)
    end
    Threads.@threads for k ∈ eachindex(indices)
        @inbounds ix, iy = indices[k]
        @inbounds x, y = xs[ix], ys[iy]
        r = _compute_radius(x, y)
        θ = mod2pi(atan(y, x))
        R = r * Δr + 1
        Θ = θ * Δθ + 1
        if R ∈ 1..nr && Θ ∈ 1..nθ
            @inbounds mc[iy,ix] = interp(Θ, R)
        end
    end
    mc
end


@inline function polar2cart(
    mp::AbstractMatrix{T};
    rows::Optional{Integer} = nothing,
    cols::Optional{Integer} = nothing,
    ν::Real = 1,
    kwargs...
) where {T <: Real,Interp <: AbstractInterp2DOrNone}
    rows = isnothing(rows) ? maybe(2 * size(mp, 2), cols) : rows
    cols = maybe(rows, cols)
    @assert ν > 0
    sθ, cθ = sincos(atan(rows, cols))
    x₀, y₀ = cθ / ν, sθ / ν
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
