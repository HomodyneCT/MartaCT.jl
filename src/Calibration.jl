module Calibration

Base.Experimental.@optlevel 3

export calibrate_image, calibrate_tomogram

using IntervalSets
using ..Monads
using ..CTImages
using ..TestImages

import Statistics; const Stat = Statistics
import ..CalibrationBase: calibrate_image, calibrate_tomogram


function calibration_helper(
    image::AbstractMatrix,
    min_pos::Union{
        Integer,
        UnitRange{<:Integer},
        CartesianIndex,
        CartesianIndices,
    },
    max_pos::Union{
        Integer,
        UnitRange{<:Integer},
        CartesianIndex,
        CartesianIndices,
    },
)
    ClosedInterval(Stat.mean.((image[min_pos], image[max_pos]))...)
end


function calibration_helper(
    image::AbstractMatrix,
    min_pos::Tuple{<:Any,<:Any},
    max_pos::Tuple{<:Any,<:Any},
)
    ClosedInterval(Stat.mean.((image[min_pos...], image[max_pos...]))...)
end


function calibration_helper(
    image::AbstractMatrix,
    min_pos::AbstractArray{C},
    max_pos::AbstractArray{C},
) where {C<:Base.AbstractCartesianIndex}
    ClosedInterval(Stat.mean.((image[min_pos], image[max_pos]))...)
end


function calibration_parameters(
    image::AbstractMatrix;
    min_pos,
    max_pos,
    interval::ClosedInterval = 0..1,
)
    cmin, cmax = calibration_helper(image, min_pos, max_pos) |> endpoints

    a, b = endpoints(interval)
    Δ = cmax - cmin
    slope = (b - a) / Δ
    intercept = (cmax * a - cmin * b) / Δ

    slope, intercept
end


"""
    calibrate_image(image::AbstractMatrix{T}; min_pos, max_pos, interval=0..1, window=nothing) where {T<:Real}

Perform calibration of `image` with reference values.
"""
function calibrate_image(
    image::AbstractMatrix{T};
    min_pos,
    max_pos,
    interval::ClosedInterval = 0..1,
    window::Optional{ClosedInterval{U}} = nothing,
) where {T<:Real,U<:Real}
    calibration = calibration_helper(image, min_pos, max_pos)
    rescale(image; interval, calibration, window)
end


function calibrate_image(img::AbstractCTImage; kwargs...)
    mmap(calibrate_image(; kwargs...), img)
end


"""
    calibrate_image(image::AbstractMatrix{T}, imp::ImageParams; interval=nothing, window=nothing) where {T<:Real}

Perform calibration of `image` using image parameters as reference.
"""
function calibrate_image(
    image::AbstractMatrix{T},
    imp::ImageParams;
    interval::Optional{ClosedInterval{U}} = nothing,
    window::Optional{ClosedInterval{W}} = nothing,
) where {T<:Real,U<:Real,W<:Real}
    min_pos = background_position(imp)
    max_pos = circle_position(imp)
    if isnothing(interval)
        interval = imp.background..imp.calibration_value
    end
    calibrate_image(
        image;
        interval,
        min_pos,
        max_pos,
        window,
    )
end


function calibrate_image(img::AbstractCTImage, imp::ImageParams; kwargs...)
    mmap(calibrate_image(imp; kwargs...), img)
end


"""
    calibrate_tomogram(image::AbstractMatrix{T}, imp::ImageParams; interval=nothing, window=nothing) where {T<:Real}

Perform calibration of reconstructed `image` with image parameters as reference.
"""
function calibrate_tomogram(
    image::AbstractMatrix{T},
    imp::ImageParams;
    interval::Optional{ClosedInterval{U}} = nothing,
    window::Optional{ClosedInterval{W}} = nothing,
) where {T<:Real,U<:Real,W<:Real}
    calibrate_image(image, imp; interval, window)
end


function calibrate_tomogram(img::CTTomogram, imp::ImageParams; kwargs...)
    mmap(calibrate_tomogram(imp; kwargs...), img)
end


"""
    calibrate_tomogram(image::AbstractMatrix{T}; min_pos, max_pos, interval=0..1, window=nothing) where {T<:Real}

Perform calibration of reconstructed `image` with reference values.
"""
function calibrate_tomogram(
    image::AbstractMatrix{T};
    min_pos,
    max_pos,
    interval::ClosedInterval = 0..1,
    window::Optional{ClosedInterval{U}} = nothing,
) where {T<:Real,U<:Real}
    calibrate_image(
        image;
        min_pos,
        max_pos,
        interval,
        window,
    )
end


function calibrate_tomogram(img::CTTomogram; kwargs...)
    mmap(calibrate_tomogram(; kwargs...), img)
end

end # module
