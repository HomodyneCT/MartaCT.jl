module CTImages

export rescale!, rescale, rotate
export sample_sinog, polar2cart
export AbstractCTImage, CTImage, CTSinogram, CTTomogram, CTImageOrTomog
export ctimage, ctsinogram, cttomogram

using IntervalSets
using ..Monads, ..Applicative, ..Interpolation
using ..Geometry: AbstractParallelBeamGeometry, AbstractFanBeamGeometry,
    ParallelBeamGeometry
import ..Monads: mbind
import ..Marta: _atype
import Base: size, map, convert


function ctimage end
function ctsinogram end
function cttomogram end


abstract type AbstractCTImage{M<:AbstractArray{<:Number}} <: Monad end

size(img::AbstractCTImage) = mbind(size, img)
size(img::AbstractCTImage, dim::Integer) = mbind(x -> size(x, dim), img)
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
    rotate(mat::AbstractMatrix{T}, α; <keyword arguments>) where {T <: Real}

Rotate matrix `mat` about the center of angle `α` given in degrees.
If `rows` and `cols` are not given, the rotated matrix has the same
dimensions of the original matrix.

# Arguments
- `mat`: matrix to rotate.
- `α`: angle in degrees.
- `rows=nothing`: number of rows of the rotated matrix.
- `cols=nothing`: number of columns of the rotated matrix.
- `interpolation`: interpolation strategy. By default is
    `LinearInterpolation`.
"""
function rotate(
    mat::AbstractMatrix{T},
    α;
    rows = nothing,
    cols = nothing,
    interpolation = mat -> LinearInterpolation(axes(mat), mat),
) where {T<:Real}
    orows, ocols = size(mat)
    if isnothing(rows)
        rows = orows
    end
    if isnothing(cols)
        cols = ocols
    end
    ϕ = deg2rad(α)
    cϕ = cos(ϕ)
    sϕ = sin(ϕ)
    rmat = zeros(T, rows, cols)
    interp = interpolation(mat)
    x₀::T = T(cols + 1) / T(2)
    y₀::T = T(rows + 1) / T(2)
    x′₀::T = T(ocols + 1) / T(2)
    y′₀::T = T(orows + 1) / T(2)
    for ix in 1:rows, iy in 1:cols
        x::T = T(ix) - x₀
        y::T = T(iy) - y₀
        x′::T = x * cϕ + y * sϕ + x′₀
        y′::T = y * cϕ - x * sϕ + y′₀
        if 1 <= x′ <= ocols && 1 <= y′ <= orows
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
    x::AbstractArray{T};
    interval::ClosedInterval = zero(T)..one(T),
    calibration::Optional{ClosedInterval{U}} = nothing,
    window::Optional{ClosedInterval{W}} = nothing,
) where {T<:Number,U<:Number,W<:Number}
    if isnothing(calibration)
        calibration = ClosedInterval(extrema(x)...)
    end

    a, b = T.(endpoints(interval))
    m, M = T.(endpoints(calibration))

    @assert m != M "Calibration values m, M = $(calibration) should be m!=M"

    f = (b - a) / (M - m)

    if m != zero(T)
        x .-= m
    end

    x .*= f

    if a != zero(T)
        x .+= a
    end

    if isnothing(window)
        return x
    end

    a, b = T.(endpoints(window))

    for i in eachindex(x)
        val = x[i]
        if val < a
            x[i] = a
        elseif val > b
            x[i] = b
        end
    end

    x
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
) where {T<:Number,U<:Number}
    res = @. image * slope + intercept

    if isnothing(window)
        return res
    end

    a, b = T.(endpoints(window))

    map(res) do x
        if x < a
            return a
        elseif x > b
            return b
        else
            return x
        end
    end
end



"""
    rescale!(x::AbstractArray, slope::Number, intercept::Number)

In place linear rescaling of `x`.

See also: [`rescale`](@ref)
"""
function rescale!(
    x::AbstractArray{T},
    slope::Number,
    intercept::Number;
    window::Optional{ClosedInterval{U}} = nothing,
) where {T<:Number,U<:Number}
    @. x = slope * x + intercept

    if isnothing(window)
        return x
    end

    a, b = T.(endpoints(window))

    for i in eachindex(x)
        val = x[i]
        if val < a
            x[i] = a
        elseif val > b
            x[i] = b
        end
    end

    x
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
    calibration  = nothing,
    window = nothing,
)
    rescale!(deepcopy(x); interval, calibration, window)
end



function sample_sinog(sinog::AbstractMatrix{T}, n) where {T<:Real}
    nd, nϕ = size(sinog)
    cumsum_sinog = rescale!(cumsum(sinog, dims = 1))
    low_sample = zeros(eltype(sinog), nd, nϕ)
    Threads.@threads for iϕ in 1:nϕ
        for _ in 1:n
            u = rand()
            for x′ in nd:-1:1
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
    background::Real = zero(T),
) where {T <: Real, Interp <: Union{Function,AbstractBilinearInterpolation}}
    nϕ, nr = size(mp)
    if isnothing(rows)
        rows = maybe(2nr, cols)
    end
    cols = maybe(rows, cols)
    mc = fill(T(background), rows, cols)

    x₀::T = T(cols - 1) / T(2)
    y₀::T = T(rows - 1) / T(2)
    Δr::T = T(1) / T(nr - 1)
    Δϕ::T = T(2π) / T(nϕ - 1)

    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(mp)

    compute_radius = (x, y) -> begin
        x == 0 && return abs(y)
        y == 0 && return abs(x)
        sqrt(x^2 + y^2)
    end

    for ix in 1:cols, iy in 1:rows
        x::T = T(ix - 1) / x₀ - T(1)
        y::T = T(iy - 1) / y₀ - T(1)
        r::T = compute_radius(x, y)
        ϕ::T = mod2pi(atan(y, x))

        ir::T = r / Δr + T(1)
        iϕ::T = ϕ / Δϕ + T(1)

        if 1 <= ir <= nr && 1 <= iϕ <= nϕ
            mc[iy, ix] = interp(iϕ, ir)
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

polar2cart(image, geometry::AbstractFanBeamGeometry; kwargs...) =
    polar2cart(image, ParallelBeamGeometry(geometry); kwargs...)

end # module
