module CTImages

export rescale!, rescale, rotate
export polar2cart
export AbstractCTImage, CTImage, CTSinogram, CTTomogram, CTImageOrTomog
export ctimage, ctsinogram, cttomogram

using IntervalSets, SimpleTraits
using ..Utils: linspace
using ..Applicative, ..Monads, ..Interpolation
using ..Geometry:
    AbstractParallelBeamGeometry, AbstractFanBeamGeometry, ParallelBeamGeometry
import ..Monads: mjoin
import ..Marta: _atype
import Base: size, convert, similar, eltype, getindex, setindex!


function ctimage end
function ctsinogram end
function cttomogram end


abstract type AbstractCTImage{M<:AbstractMatrix} <: AbstractMatrix{eltype(M)} end

@traitimpl Monad{AbstractCTImage}

mjoin(img::AbstractCTImage) = img.data

_atype(img::M) where {M <: AbstractCTImage} = _atype(M)
eltype(::Type{M}) where {M<:AbstractCTImage{A}} where A = eltype(A)
size(img::AbstractCTImage) = mbind(size, img)
size(img::AbstractCTImage, dim::Int) = size(mjoin(img), dim)
getindex(img::AbstractCTImage, i::Int) = getindex(mjoin(img), i)
getindex(img::AbstractCTImage, I::Vararg{Int,N}) where N =
    getindex(mjoin(img), I...)
setindex!(img::AbstractCTImage, v, i::Int) = setindex!(mjoin(img), v, i)
setindex!(img::AbstractCTImage, v, I::Vararg{Int,N}) where N =
    setindex!(mjoin(img), v, I...)


const ctnames = (image = :CTImage, sinog = :CTSinogram, tomog = :CTTomogram)
const ctfn = (image = :ctimage, sinog = :ctsinogram, tomog = :cttomogram)


for nm in ctnames
    @eval begin
        struct $nm{M} <: AbstractCTImage{M}
            data::M
        end
        _atype(::Type{$nm{M}}) where M = M
        $nm{M}(img::$nm) where {M} = mreturn($nm{M}, convert(M, mjoin(img)))
        convert(::Type{M}, img::$nm) where {M<:$nm} = mreturn(M, img)
        similar(::Type{M}, args...) where {M<:$nm} =
            mreturn(M, similar(_atype(T), args...))
    end
end


const CTImageOrTomog{M} = Union{CTImage{M},CTTomogram{M}}


"""
    rotate(mat::AbstractMatrix{T}, α::Real; <keyword arguments>) where {T <: Real}

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
    rows::Optional{<:Integer} = nothing,
    cols::Optional{<:Integer} = nothing,
    background::Optional{<:Real} = nothing,
    interpolation::Optional{<:Interp} = nothing,
) where {T <: Real,Interp <: AbstractInterp2DOrNone}
    orows, ocols = size(mat)
    rows = maybe(orows, rows)
    cols = maybe(ocols, cols)
    sϕ, cϕ = sincos(deg2rad(α))
    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(mat)
    x₀::T = T(cols + 1) / 2
    y₀::T = T(rows + 1) / 2
    x′₀::T = T(ocols + 1) / 2
    y′₀::T = T(orows + 1) / 2
    z::T = maybe(zero(T), background)
    rmat = similar(mat, rows, cols)
    fill!(rmat, z)
    for ix ∈ 1:rows, iy ∈ 1:cols
        x::T = T(ix) - x₀
        y::T = T(iy) - y₀
        x′::T = x * cϕ + y * sϕ + x′₀
        y′::T = y * cϕ - x * sϕ + y′₀
        if x′ ∈ 1..ocols && y′ ∈ 1..orows
            rmat[iy, ix] = interp(y′, x′)
        end
    end
    rmat
end


"""
    rescale!(x::AbstractArray{T,N}; interval=nothing, calibration=nothing, window=nothing) where {T,N}

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
    map(res) do x
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
    map(img) do x
        x < a && return a
        x > b && return b
        return x
    end
    img
end



"""
    rescale(x::AbstractArray{T,N}; interval=nothing, calibration=nothing, window=nothing) where {T,N}

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


for nm ∈ (:rescale, :rescale!)
    @eval begin
        $nm(img::AbstractCTImage, args...) = mmap(img) do x
            $nm(x, args...)
        end
        $nm(img::AbstractCTImage; kwargs...) = mmap(img) do x
            $nm(x; kwargs...)
        end
    end
end


function polar2cart(
    mp::AbstractMatrix{T};
    rows::Optional{<:Integer} = nothing,
    cols::Optional{<:Integer} = nothing,
    background::Optional{<:Real} = nothing,
    ν::Real = 1,
    transposed::Bool = false,
    interpolation::Optional{Interp} = nothing,
) where {T <: Real,Interp <: AbstractInterp2DOrNone}
    mp = transposed ? permutedims(mp) : mp
    nϕ, nr = size(mp)
    rows = isnothing(rows) ? maybe(2nr, cols) : rows
    cols = maybe(rows, cols)
    ν = T(ν)
    X::T, Y::T = cols, rows
    y₀, x₀ = sincos(atan(Y, X))
    Δr::T = (nr - 1) / ν
    Δϕ::T = (nϕ - 1) / 2π
    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(mp)
    compute_radius = (x::T, y::T) -> begin
        x == 0 && return abs(y)
        y == 0 && return abs(x)
        return √(x^2 + y^2)
    end
    xs = linspace(-x₀, x₀, cols)
    ys = linspace(-y₀, y₀, rows)
    indices = Vector{NTuple{2,Int}}(undef, rows * cols)
    @inbounds for k ∈ eachindex(indices)
        indices[k] = (k - 1) ÷ rows + 1, (k - 1) % rows + 1
    end
    z::T = maybe(zero(T), background)
    mc = similar(mp, rows, cols)
    fill!(mc, z)
    Threads.@threads for k ∈ eachindex(indices)
        @inbounds ix, iy = indices[k]
        @inbounds x, y = xs[ix], ys[iy]
        r = compute_radius(x, y)
        ϕ = mod2pi(atan(y, x))
        R = r * Δr + 1
        Θ = ϕ * Δϕ + 1
        if R ∈ 1..nr && Θ ∈ 1..nϕ
            @inbounds mc[iy,ix] = interp(Θ, R)
        end
    end
    mc
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
