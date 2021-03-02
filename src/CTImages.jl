module CTImages

export rescale!, rescale, rotate
export sample_sinog, polar2cart, cart2polar
export AbstractCTImage, CTImage, CTSinogram, CTTomogram, CTImageOrTomog
export ctimage, ctsinogram, cttomogram

using IntervalSets
using ..Marta:linspace
using ..Monads, ..Applicative, ..Interpolation
using ..Geometry: AbstractParallelBeamGeometry, AbstractFanBeamGeometry,
    ParallelBeamGeometry
import ..Monads:mbind
import ..Marta:_atype
import Base: size, map, convert, length


function ctimage end
function ctsinogram end
function cttomogram end


abstract type AbstractCTImage{M <: AbstractArray{<:Number}} <: Monad end

size(img::AbstractCTImage) = mbind(size, img)
size(img::AbstractCTImage, dim::Integer) = mbind(x -> size(x, dim), img)
length(img::AbstractCTImage) = mbind(length, img)
map(f::Function, img::AbstractCTImage) = mmap(fmap(f), img)

_atype(::Type{<:AbstractCTImage{M}}) where {M} = M
_atype(img::AbstractCTImage) = _atype(typeof(img))


const ctnames = (image = :CTImage, sinog = :CTSinogram, tomog = :CTTomogram)
const ctfn = (image = :ctimage, sinog = :ctsinogram, tomog = :cttomogram)


for nm in ctnames
    @eval begin
        struct $nm{M} <: AbstractCTImage{M}
        data::M
    end
        $nm{M}(img::$nm) where {M} = $nm{M}(convert(M, img.data))
        mbind(f::Union{Function,Type}, m::$nm) = f(m.data)
        convert(::Type{$nm{M}}, img::$nm) where {M} = $nm{M}(convert(M, img.data))
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
    interpolation::Optional{<:Interp} = nothing,
) where {T <: Real,Interp <: Union{Function,AbstractInterp2DOrNone}}
    orows, ocols = size(mat)
    rows = maybe(orows, rows)
    cols = maybe(ocols, cols)
    sϕ, cϕ = sincos(deg2rad(α))
    rmat = zeros(T, rows, cols)
    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(mat)
    x₀::T = T(cols + 1) / 2
    y₀::T = T(rows + 1) / 2
    x′₀::T = T(ocols + 1) / 2
    y′₀::T = T(orows + 1) / 2
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
    @assert(
        m != M,
        "Cannot calibrate image as calibration values " *
        "m, M = $(calibration) are equal"
    )
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


function sample_sinog(sinog::AbstractMatrix{T}, n) where {T <: Real}
    nd, nϕ = size(sinog)
    cumsum_sinog = rescale!(cumsum(sinog, dims = 1))
    low_sample = zeros(eltype(sinog), nd, nϕ)
    Threads.@threads for iϕ ∈ 1:nϕ
            foreach(1:n) do _
            u = rand()
            for x′ ∈ nd:-1:1
                if cumsum_sinog[x′, iϕ] < u
                    low_sample[x′, iϕ] += 1
                    break
                end
            end
        end
    end
    rescale!(low_sample)
end


sample_sinog(n) = x -> sample_sinog(x, n)
sample_sinog(sinog::CTSinogram, n) = CTSinogram(sinog ↣ sample_sinog(n))


function polar2cart(
    mp::AbstractMatrix{T};
    rows::Optional{<:Integer} = nothing,
    cols::Optional{<:Integer} = nothing,
    interpolation::Optional{Interp} = nothing,
    background::Optional{<:Real} = nothing,
    ν::Real = 1,
    transposed::Bool = false,
) where {
    T <: Real,
    Interp <: Union{Function,AbstractInterp2DOrNone},
}
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
    compute_radius = (x, y) -> begin
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
    mc = fill(z, rows, cols)
    Threads.@threads for k ∈ eachindex(indices)
        @inbounds ix, iy = indices[k]
        @inbounds x, y = xs[ix], ys[iy]
        r::T = compute_radius(x, y)
        ϕ::T = mod2pi(atan(y, x))
        R::T = r * Δr + 1
        Θ::T = ϕ * Δϕ + 1
        if R ∈ 1..nr && Θ ∈ 1..nϕ
            @inbounds mc[iy,ix] = interp(Θ, R)
         end
    end
    mc
end


polar2cart(; kwargs...) = x -> polar2cart(x; kwargs...)
polar2cart(g::AbstractParallelBeamGeometry; kwargs...) = x -> polar2cart(x, g; kwargs...)

polar2cart(image::AbstractMatrix, geometry::AbstractParallelBeamGeometry; kwargs...) =
    polar2cart(image; geometry.rows, geometry.cols, kwargs...)

polar2cart(image::CTImageOrTomog, geometry::AbstractParallelBeamGeometry; kwargs...) =
    typeof(image)(image ↣ polar2cart(geometry; kwargs...))

polar2cart(image::CTImageOrTomog; kwargs...) =
    typeof(image)(image ↣ polar2cart(; kwargs...))

polar2cart(image, geometry::AbstractFanBeamGeometry; kwargs...) =
    polar2cart(image, ParallelBeamGeometry(geometry); kwargs...)

end # module
