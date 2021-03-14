module TestImages

export indicator_matrix, indicator_matrix_2
export AbstractImageParams, ImageParams, CircleParams
export circle_position, background_position
export gray_scale_indices
export circle_image, gray_scale_image
export square_image, circle_polar_image
export combined_images_size, combine_images
export create_image, plateau_length
export AbstractTestImage, AbstractGrayScale
export GrayScaleLine, GrayScalePyramid, CircleImage, WhiteRect


using IntervalSets
using ..Monads
using ..CTImages: CTImage, polar2cart
import ..CTImages: rescale, rescale!
using ..Interpolation: NoInterpolation
import ..Geometry: ParallelBeamGeometry, FanBeamGeometry, AbstractParallelBeamGeometry
import ..Geometry
using ..AbstractAlgorithms: AbstractProjectionAlgorithm
import ..Marta: datatype, yaml_repr, struct2dict
import ..Info: CTInfo
import ..AbstractAlgorithms: project_image
import ..CalibrationBase: calibrate_image, calibrate_tomogram, calibration_data
import Base: getproperty, size, show


abstract type AbstractImageParams end

size(p::AbstractImageParams) = p.rows, p.cols
function size(p::AbstractImageParams, dim::Integer)
    dim == 1 && return p.rows
    dim == 2 && return p.cols
    return 1
end

yaml_repr(imp::AbstractImageParams) = struct2dict(imp)
CTInfo(imp::AbstractImageParams) = CTInfo(pairs(struct2dict(imp))...)


function indicator_matrix(::Type{T}, sq::Int, pad::Int = 0) where {T<:Real}
    if pad == 0
        pad = sq + 1
    end
    rows = 2 * (sq + pad)
    cols = 2 * (sq + pad)
    matrix = zeros(T, rows, cols)
    matrix[pad+sq+1:pad+2sq, pad+1:pad+sq] .= 1.0
    matrix[pad+1:pad+sq, pad+sq+1:pad+2sq] .= 1.0
    return matrix
end

indicator_matrix(sq::Int) = indicator_matrix(Float32, sq)


function indicator_matrix_2(
    ::Type{T},
    sq1::Int,
    sq2::Int,
    pad::Int = 0,
) where {T<:Real}
    sqmin, sqmax = minmax(sq1, sq2)
    if pad == 0
        pad = sqmin + 1
    end
    rows = 2 * (sqmax + pad)
    cols = 2 * (sqmax + pad)
    matrix = zeros(T, rows, cols)
    matrix[pad+1:pad+sqmin, pad+1:pad+sqmin] .= 1.0
    matrix[pad+sqmax+1:pad+2sqmax, pad+sqmax+1:pad+2sqmax] .= 1.0
    return matrix
end

indicator_matrix_2(sq1::Int, sq2::Int, pad::Int = 0) =
    indicator_matrix_2(Float32, sq1, sq2, pad)


const _default_gray_scale_params = Dict(
    :swidth => 200, # Gray scale width
    :sheight => 40, # Gray scale height
    :pad => 30,     # Padding around gray scale and calibration circle
    :dist => 10,    # Distance between gray scale and calibration circle
    :radius => 15,  # Calibration circle radius
)


"""
    struct ImageParams{T<:Real}

Hold information to construct the input test image.
"""
struct ImageParams{T<:Real} <: AbstractImageParams
    rows::Int # Test image rows
    cols::Int # Test image cols
    gray_scale::ClosedInterval{T} # Gray scale interval.
    calibration_value::T # Value for algorithm calibration.
    background::T
    kw::Dict{Symbol,Any}

    function ImageParams{T}(
        rows::Optional{<:Integer} = nothing,
        cols::Optional{<:Integer} = nothing,
        gray_scale::ClosedInterval = -1000..1000,
        calibration_value::Optional{<:Real} = nothing,
        background::Optional{<:Real} = nothing;
        kwargs...
    ) where {T<:Real}
        kw = isempty(kwargs) ? _default_gray_scale_params : Dict(kwargs)
        srows, scols = combined_images_size(; kw...)
        rows, cols = maybe(srows, rows), maybe(scols, cols)
        @assert(
            rows >= srows && cols >= scols,
            "Requested image size is too small: " *
            "rows($rows) ≥ $srows && columns($cols) ≥ $scols"
        )
        background = maybe(leftendpoint(gray_scale), background)
        calibration_value = isnothing(calibration_value) ?
            mean(gray_scale) : calibration_value

        if calibration_value ∉ background..rightendpoint(gray_scale)
            midpoint = mean(gray_scale)
            @warn(
                "Calibration value should be in the interval of the gray scale ($calibration_value ∉ $gray_scale), taking midpoint: $midpoint"
            )
            calibration_value = midpoint
        end

        new(
            rows,
            cols,
            gray_scale,
            calibration_value,
            background,
            kw,
        )
    end
end


datatype(imp::ImageParams{T}) where T = T


function getproperty(p::ImageParams, s::Symbol)
    s ∈ fieldnames(ImageParams) && return getfield(p, s)
    s === :width && return p.cols
    s === :height && return p.rows
    p.kw[s]
end


"""
    ImageParams([T=Float32]; <keyword arguments>) where {T<:Real}

Construct ImageParams object.

# Arguments:
- `width=200`: width of the rectangle gray scale image.
- `height=40`: height of the rectangle gray scale image.
- `pad=30`: padding around image.
- `dist=10`: distance between images.
- `radius=20`: radius of calibration circle.
- `gray_scale=0..1`: interval of gray scale values.
- `calibration_value=nothing`: value of the calibration circle. If not specified, is `gray_scale[2]`.
"""
function ImageParams(
    ::Type{T} = Float32;
    rows::Optional{<:Integer} = nothing,
    cols::Optional{<:Integer} = nothing,
    width::Optional{<:Integer} = nothing,
    height::Optional{<:Integer} = nothing,
    gray_scale::ClosedInterval = -1000..1000,
    calibration_value::Optional{<:Real} = nothing,
    background::Optional{<:Real} = nothing,
    kwargs...
) where {T<:Real}
    rows = maybe(rows, height)
    cols = maybe(cols, width)
    ImageParams{T}(
        rows,
        cols,
        gray_scale,
        calibration_value,
        background;
        kwargs...
    )
end


show(io::IO, p::AbstractImageParams) = print(
    io,
    """Image Params:
- size (W×H): $(p.width) × $(p.height)
- gray scale: $(p.gray_scale)
- calibration value: $(p.calibration_value)
- background value: $(p.background)""",
)


struct CircleParams{T <: Real} <: AbstractImageParams
    radius::Int # Calibration circle radius
    rows::Int # Test image rows
    cols::Int # Test image cols
    gray_scale::ClosedInterval{T} # Gray scale interval.
    calibration_value::T # Value for algorithm calibration.
    background::T

    function CircleParams{T}(
        radius::Integer,
        rows::Integer,
        cols::Integer,
        gray_scale::ClosedInterval = -1000..1000,
        calibration_value::Optional{<:Real} = nothing,
        background::Optional{<:Real} = nothing,
    ) where {T<:Real}
        background = maybe(leftendpoint(gray_scale), background)
        calibration_value = isnothing(calibration_value) ?
            mean(gray_scale) : calibration_value

        if calibration_value ∉ background..rightendpoint(gray_scale)
            midpoint = mean(gray_scale)
            @warn "Calibration value should be in the interval of the gray scale ($calibration_value ∉ $gray_scale), taking midpoint: $midpoint"
            calibration_value = midpoint
        end

        new(
            radius,
            rows,
            cols,
            gray_scale,
            calibration_value,
            background,
        )
    end
end


datatype(imp::CircleParams{T}) where T = T

function getproperty(p::CircleParams, s::Symbol)
    s ≡ :width && return p.cols
    s ≡ :height && return p.rows
    getfield(p, s)
end


function CircleParams(
    ::Type{T} = Float32;
    radius::Optional{<:Integer} = nothing,
    width::Optional{<:Integer} = nothing,
    height::Optional{<:Integer} = nothing,
    rows::Optional{<:Integer} = nothing,
    cols::Optional{<:Integer} = nothing,
    gray_scale::ClosedInterval = -1000..1000,
    calibration_value::Optional{<:Real} = nothing,
    background::Optional{<:Real} = nothing,
) where {T <: Real}
    factor = 31 / 50 # magic number!
    circ_zoom = 1.3 / 17
    #=
        If width or height are present, use them, otherwise use rows and cols.
        Finally, if rows or cols are not present, compute defaults.
    =#
    rows, cols = maybe(rows, height), maybe(cols, width)
    width = maybe(500, cols)
    height = isnothing(rows) ? round(Int, width * factor) : rows
    radius = isnothing(radius) ? round(Int, width * circ_zoom) : radius
    rows, cols = height, width
    CircleParams{T}(radius, rows, cols, gray_scale, calibration_value, background)
end


"""
    circle_position(imp::AbstractImageParams)

Return circle position inside image given image parameters.
"""
function circle_position end

function circle_position(imp::ImageParams)
    srows, _ = combined_images_size(imp)
    rows, cols = size(imp)
    r = (rows - srows) ÷ 2 + imp.pad + imp.radius + 1
    c = cols ÷ 2
    r, c
end
circle_position(imp::CircleParams) = imp.rows ÷ 2, imp.cols ÷ 2


"""
    background_position(imp::AbstractImageParams)

Return suitable position to calibrate background.
"""
function background_position end

function background_position(imp::ImageParams)
    cr, cc = circle_position(imp)
    cr, cc ÷ 2
end

function background_position(imp::CircleParams)
    cr, cc = circle_position(imp)
    max(1, cr - imp.radius), max(1, cc - imp.radius)
end


"""
    gray_scale_indices(imp::ImageParams)

Return a tuple `(row, column_range)` where `row` is the row index of the scale and `column_range` is the range of columns.
"""
function gray_scale_indices(imp::ImageParams)
    rows, cols = size(imp)
    srows, scols = combined_images_size(imp)
    min_row = (rows - srows) ÷ 2 + imp.pad + 2 * imp.radius + imp.dist + 1
    max_row = min_row + imp.sheight - 1
    min_col = (cols - scols) ÷ 2 + imp.pad + 1
    max_col = min_col + imp.swidth - 1
    min_row:max_row, min_col:max_col
end


"""
    circle_image([T=Float32]; radius=20, value=1) where {T<:Real}

Create a square image with a circle of given value.

# Examples
```julia-repl
julia> circle_image(30, 0.8)
30×30 Array{Float32,2}:
[...]
```

See also: [`gray_scale_image`](@ref), [`combine_images`](@ref)
"""
function circle_image(
    ::Type{T} = Float32;
    radius::Integer = 20,
    calibration_value::Real = zero(T),
    background::Real = -1000,
    rows::Optional{<:Integer} = nothing,
    cols::Optional{<:Integer} = nothing,
    swidth::Optional{<:Integer} = nothing,
    sheight::Optional{<:Integer} = nothing,
    ν::Optional{<:Real} = nothing,
    kwargs...,
) where {T<:Real}
    rows = isnothing(sheight) ? maybe(2radius, rows) : sheight
    cols = isnothing(swidth) ? maybe(rows, cols) : swidth
    nr = min(rows, cols) ÷ 2
    ν = isnothing(ν) ? min(rows,cols) / hypot(rows, cols) : ν
    polar2cart(
        circle_polar_image(T; nr, radius, calibration_value, background);
        rows, cols, background, ν, kwargs...)
end


"""
    circle_image(imp::AbstractImageParams)

Create a square image with a circle of given value from parameters `imp`.
"""
function circle_image end

function circle_image(imp::ImageParams; kwargs...)
    circle_image(
        datatype(imp);
        imp.radius,
        imp.calibration_value,
        imp.background,
        rows = 2imp.radius,
        kwargs...,
    )
end

function circle_image(imp::CircleParams; kwargs...)
    circle_image(
        datatype(imp);
        imp.radius,
        imp.calibration_value,
        imp.background,
        imp.rows,
        imp.cols,
        kwargs...,
    )
end


"""
    gray_scale_image(
        [T=Float32];
        rows=40,
        cols=200,
        swidth=nothing,
        sheight=nothing,
        gray_scale=-1000..1000,
    ) where {T <: Real}

Create an image with a gray scale rectangle with given scale gray_scale.

# Examples
```julia-repl
julia> gray_scale_image(swidth=80, sheight=40, gray_scale=-1000..1000)
40×80 Array{Float32,2}:
[...]
```

See also: [`circle_image`](@ref), [`combine_images`](@ref)
"""
function gray_scale_image(
    ::Type{T} = Float32;
    rows::Integer = 200,
    cols::Integer = 40,
    swidth::Optional{<:Integer} = nothing,
    sheight::Optional{<:Integer} = nothing,
    gray_scale::ClosedInterval = -1000..1000,
    background::Optional{<:Real} = nothing,
) where {T<:Real}
    rows = maybe(rows, sheight)
    cols = maybe(cols, swidth)
    minv, maxv = endpoints(gray_scale) .|> T
    background = maybe(infimum(gray_scale), background)
    val_range = minv == maxv ? minv : range(minv, maxv, length = cols)
    image = fill(T(background), rows, cols)
    image[1:rows,1:cols] .= val_range'
    image
end


"""
pyramid_gray_scale_image([T=Float32]; swidth=200, sheight=40, gray_scale=-1000..1000, plateau = 0) where {T <: Real}

Create an image with a pyramid gray scale rectangle with given scale gray_scale.

# Examples
```julia-repl
julia> pyramid_gray_scale_image(swidth=80, sheight=40, gray_scale=-1000..1000)
40×80 Array{Float32,2}:
[...]
```

See also: [`gray_scale_image`](@ref), [`circle_image`](@ref), [`combine_images`](@ref)
"""
function pyramid_gray_scale_image(
    ::Type{T} = Float32;
    swidth::Integer = 200,
    sheight::Integer = 40,
    gray_scale::ClosedInterval = -1000..1000,
    background::Optional{<:Real} = nothing,
    plateau::Real = 0,
) where {T<:Real}
    @assert 0 ≤ plateau ≤ 1 "Plateau is a fraction of the gray scale width!"
    minv, maxv = endpoints(gray_scale) .|> T
    val_range = fill(maxv, swidth)
    plateau_len = round(Int, swidth * plateau)
    len = (swidth - plateau_len) ÷ 2
    if len ≠ 0
        step = (maxv - minv) / (swidth - plateau_len - len - 1)
        half_scale = range(minv; step, length = len)
        val_range[1:len] .= half_scale
        val_range[end:-1:(swidth-len+1)] .= half_scale # Odd case included.
    end
    background = maybe(minv, background)
    image = fill(background, sheight, swidth)
    image .= val_range'
    image
end


"""
    gray_scale_image(imp::ImageParams)

Create an image with a gray scale rectangle with given scale gray_scale
from parameters `imp`.
"""
function gray_scale_image(imp::ImageParams{T}) where {T<:Real}
    gray_scale_image(
        T; imp.swidth, imp.sheight, imp.gray_scale, imp.background)
end


"""
    pyramid_gray_scale_image(imp::ImageParams)

Create an image with a pyramid gray scale rectangle with given scale gray_scale
from parameters `imp`.

Create an image with a pyramid gray scale rectangle with given scale gray_scale
from parameters `imp`.
"""
function pyramid_gray_scale_image(
    imp::ImageParams{T};
    plateau::Real = 0,
) where {T<:Real}
    pyramid_gray_scale_image(
        T;
        imp.swidth,
        imp.sheight,
        imp.gray_scale,
        imp.background,
        plateau,
    )
end


function combined_images_size(; swidth, sheight, radius, pad, dist)
    sheight + 2radius + 2pad + dist, swidth + 2pad
end

function combined_images_size(p::ImageParams)
    combined_images_size(;
        p.swidth,
        p.sheight,
        p.radius,
        p.pad,
        p.dist,
    )
end


"""
    combine_images(imp::ImageParams{T}, rect::AbstractMatrix{T}, circle::AbstractMatrix{T}) where {T <: Real}

Helper function to create a gray scale image.

Combine the `rect` image with the `circle` image for calibration.
The circle image should be smaller than the gray scale image!

See also: [`gray_scale_image`](@ref), [`circle_image`](@ref)
"""
function combine_images(
    imp::ImageParams{T},
    rect::AbstractMatrix{T},
    circle::AbstractMatrix{T},
) where {T<:Real}
    rect_rows, rect_cols = size(rect)
    circ_rows, circ_cols = size(circle)

    @assert((rect_rows, rect_cols) == (imp.sheight, imp.swidth),
        "Expected gray scale size $((imp.sheight, imp.swidth)), got $((srows, scols)) instead")

    @assert(2imp.radius == circ_rows,
        "Expected radius $(imp.radius), got $(circ_rows ÷ 2)")

    if circ_rows >= rect_rows || circ_cols >= rect_cols
        @warn(
            "Circle image should be smaller than gray scale rectangle",
            (circ_rows, circ_cols),
            (rect_rows, rect_cols)
        )
    end

    srows, scols = combined_images_size(imp)
    rows, cols = size(imp)
    @assert(
        rows >= srows && cols >= scols,
        "Requested image size too small: " *
        "rows($rows) ≥ $srows && columns($cols) ≥ $scols"
    )

    simage = fill(imp.background, srows, scols)
    image = fill(imp.background, rows, cols)

    crow_beg = imp.pad + 1
    crow_end = crow_beg + circ_rows - 1
    ccol_beg = scols ÷ 2 - imp.radius + 1
    ccol_end = ccol_beg + circ_cols - 1

    simage[crow_beg:crow_end, ccol_beg:ccol_end] .= circle

    rrow_beg = crow_end + imp.dist + 1
    rrow_end = rrow_beg + rect_rows - 1
    rcol_beg = imp.pad + 1
    rcol_end = rcol_beg + rect_cols - 1

    simage[rrow_beg:rrow_end, rcol_beg:rcol_end] .= rect

    ri = (rows - srows) ÷ 2 + 1
    rf = ri + srows - 1
    ci = (cols - scols) ÷ 2 + 1
    cf = ci + scols - 1
    image[ri:rf,ci:cf] .= simage
    image
end


"""
    create_image(par::ImageParams)

Create gray scale image.
"""
function create_image(par::ImageParams)
    circ_img = circle_image(par)
    rect_img = gray_scale_image(par)
    combine_images(par, rect_img, circ_img) |> CTImage
end


"""
    create_pyramid_image(par::ImageParams)

Create pyramid gray scale image.
"""
function create_pyramid_image(par::ImageParams; plateau::Real = 0)
    circ_img = circle_image(par)
    rect_img = pyramid_gray_scale_image(par; plateau)
    combine_images(par, rect_img, circ_img) |> CTImage
end


abstract type AbstractTestImage end
abstract type AbstractGrayScale <: AbstractTestImage end

size(gs::AbstractTestImage) = size(gs.image)
size(gs::AbstractTestImage, d::Integer) = size(gs.image, d)


struct CircleImage{T<:Real} <: AbstractTestImage
    params::CircleParams{T}
    image::CTImage{Matrix{T}}
end


function getproperty(grimg::CircleImage, s::Symbol)
    s ∈ fieldnames(CircleImage) && return getfield(grimg, s)
    getproperty(grimg.params, s)
end


function CircleImage(imp::CircleParams; kwargs...)
    image = circle_image(imp; kwargs...)
    CircleImage(imp, CTImage(image))
end


function CircleImage(
    ::Type{T} = Float32; interpolation = nothing, kwargs...) where{T <: Real}
    interpolation = maybe(NoInterpolation(), interpolation)
    CircleImage(CircleParams(T; kwargs...); interpolation)
end


function CircleImage(imp::ImageParams; kwargs...)
    CircleImage(
        datatype(imp);
        imp.radius,
        imp.rows,
        imp.cols,
        imp.gray_scale,
        imp.calibration_value,
        imp.background,
        kwargs...
    )
end


function CircleImage(pbg::AbstractParallelBeamGeometry; kwargs...)
    CircleImage(datatype(pbg); pbg.width, pbg.height, kwargs...)
end

struct GrayScaleLine{T<:Real} <: AbstractGrayScale
    params::ImageParams{T}
    image::CTImage{Matrix{T}}
end


function getproperty(grimg::GrayScaleLine, s::Symbol)
    s ∈ fieldnames(GrayScaleLine) && return getfield(grimg, s)
    getproperty(grimg.params, s)
end


function GrayScaleLine(imp::ImageParams; kwargs...)
    GrayScaleLine(imp, create_image(imp); kwargs...)
end


function GrayScaleLine(
    ::Type{T} = Float32;
    radius::Optional{<:Integer} = nothing,
    width::Optional{<:Integer} = nothing,
    height::Optional{<:Integer} = nothing,
    rows::Optional{<:Integer} = nothing,
    cols::Optional{<:Integer} = nothing,
    gray_scale::ClosedInterval = -1000..1000,
    calibration_value::Real = zero(T),
    background::Real = -1000,
    kwargs...,
) where {T<:Real}
    factor = 31 / 50 # magic number!
    #=
        If width or height are present, use them, otherwise use rows and cols.
        Finally, if rows or cols are not present, compute defaults.
    =#
    rows, cols = maybe(rows, height), maybe(cols, width)
    width = maybe(500, cols)
    height = isnothing(rows) ? round(Int, width * factor) : rows
    circ_zoom = 1
    swidth = round(Int, 17 / 25 * width)
    sheight = round(Int, 6 / 31 * height)
    pad = round(Int, 3 / 25 * width / circ_zoom)
    dist = 2pad ÷ 3
    radius = isnothing(radius) ? round(Int, width * circ_zoom / 17) : radius
    rows, cols = height, width
    imp = ImageParams(
        T;
        rows,
        cols,
        swidth,
        sheight,
        pad,
        dist,
        radius,
        gray_scale,
        calibration_value,
        background,
    )
    GrayScaleLine(imp; kwargs...)
end


function GrayScaleLine(geometry::AbstractParallelBeamGeometry; kwargs...)
    GrayScaleLine(
        datatype(geometry);
        geometry.width,
        geometry.height,
        kwargs...
    )
end


struct GrayScalePyramid{T<:Real} <: AbstractGrayScale
    params::ImageParams{T}
    plateau::T
    image::CTImage{Matrix{T}}
end


function getproperty(grimg::GrayScalePyramid, s::Symbol)
    s ∈ fieldnames(GrayScalePyramid) && return getfield(grimg, s)
    getproperty(grimg.params, s)
end


function GrayScalePyramid(
    imp::ImageParams{T}; plateau::Real = zero(T), kwargs...) where {T<:Real}
    image = create_pyramid_image(imp; plateau)
    GrayScalePyramid(imp, plateau, image)
end


function GrayScalePyramid(
    ::Type{T} = Float32;
    radius::Optional{<:Integer} = nothing,
    width::Optional{<:Integer} = nothing,
    height::Optional{<:Integer} = nothing,
    rows::Optional{<:Integer} = nothing,
    cols::Optional{<:Integer} = nothing,
    gray_scale::ClosedInterval = -1000..1000,
    calibration_value::Real = zero(T),
    background::Real = -1000,
    plateau::Real = zero(T),
    kwargs...,
) where {T<:Real}
    factor = 31 / 50 # magic number!
    #=
        If width or height are present, use them, otherwise use rows and cols.
        Finally, if rows or cols are not present, compute defaults.
    =#
    rows, cols = maybe(rows, height), maybe(cols, width)
    width = maybe(500, cols)
    height = isnothing(rows) ? round(Int, width * factor) : rows
    circ_zoom = 1
    swidth = round(Int, 17 / 25 * width)
    sheight = round(Int, 6 / 31 * height)
    pad = round(Int, 3 / 25 * width / circ_zoom)
    dist = 2pad ÷ 3
    radius = isnothing(radius) ? round(Int, width * circ_zoom / 17) : radius
    rows, cols = height, width
    imp = ImageParams(
        T;
        rows,
        cols,
        swidth,
        sheight,
        pad,
        dist,
        radius,
        gray_scale,
        calibration_value,
        background,
    )
    GrayScalePyramid(imp; plateau, kwargs...)
end


function GrayScalePyramid(geometry::AbstractParallelBeamGeometry; kwargs...)
    GrayScalePyramid(
        datatype(geometry);
        geometry.width,
        geometry.height,
        kwargs...
    )
end


"""
    circle_polar_image([T=Float32] nr, nϕ, radius; value=1) where {T<:Real}

Create a circle of radius `radius` inside a `nr×nϕ` image
in polar coordinates .

# Examples
```julia-repl
julia> circle_polar_image(30, 30, 10)
30×30 Array{Float32,2}:
[...]
```
"""
function circle_polar_image(
    ::Type{T},
    nr::Integer,
    nϕ::Integer,
    radius::Integer;
    calibration_value::Real = zero(T),
    background::Real = -1000,
) where {T<:Real}
    @assert radius <= nr "Radius should be less than number of radial points"
    mat = fill(T(background), nϕ, nr)
    mat[:, 1:radius] .= calibration_value
    mat
end
circle_polar_image(::Type{T}, nr::Integer, radius::Integer; kwargs...) where {T} =
    circle_polar_image(T, nr, nr, radius; kwargs...)
circle_polar_image(nr::Integer, nϕ::Integer, radius::Integer; kwargs...) =
    circle_polar_image(Float32, nr, nϕ, radius; kwargs...)
circle_polar_image(nr::Integer, radius::Integer; kwargs...) =
    circle_polar_image(nr, nr, radius; kwargs...)
circle_polar_image(::Type{T}, radius::Integer; kwargs...) where {T} =
    circle_polar_image(T, radius, radius; kwargs...)
circle_polar_image(radius::Integer; kwargs...) =
    circle_polar_image(radius, radius; kwargs...)
function circle_polar_image(
    ::Type{T} = Float32;
    radius::Integer = 60,
    nr::Optional{<:Integer} = nothing,
    nϕ::Optional{<:Integer} = nothing,
    kwargs...) where {T<:Real}
    nr = maybe(2radius, nr)
    nϕ = maybe(nr, nϕ)
    circle_polar_image(T, nr, nϕ, radius; kwargs...)
end


"""
    square_image([T=Float32] r, c; l=nothing) where {T<:Real}

Create a `l×l` square inside a `r×c` image.

# Examples
```julia-repl
julia> square_image(30, 30; l=10)
30×30 Array{Float32,2}:
[...]
```
"""
function square_image(
    ::Type{T},
    r::Integer,
    c::Integer;
    l::Optional{<:Integer} = nothing,
    calibration_value::Real = zero(T),
    background::Real = -1000
) where {T<:Real}
    matrix = fill(T(background), r, c)

    if isnothing(l)
        l = min(r, c)
    end

    @assert l <= min(r, c)

    rin = r ÷ 2 - l ÷ 2 + 1
    rfi = rin + l - 1
    cin = c ÷ 2 - l ÷ 2 + 1
    cfi = cin + l - 1

    matrix[rin:rfi, cin:cfi] .= calibration_value

    matrix
end


square_image(r, c; kwargs...) = square_image(Float32, r, c; kwargs...)


struct WhiteRect{T<:Real} <: AbstractGrayScale
    params::ImageParams{T}
    image::CTImage{Matrix{T}}
end


function getproperty(grimg::WhiteRect, s::Symbol)
    s ∈ fieldnames(WhiteRect) && return getfield(grimg, s)
    getproperty(grimg.params, s)
end


function WhiteRect(imp::ImageParams; kwargs...)
    grimg = GrayScaleLine(imp; kwargs...)
    WhiteRect(grimg.params, grimg.image)
end


function WhiteRect(
    ::Type{T} = Float32;
    radius::Optional{<:Integer} = nothing,
    width::Optional{<:Integer} = nothing,
    height::Optional{<:Integer} = nothing,
    rows::Optional{<:Integer} = nothing,
    cols::Optional{<:Integer} = nothing,
    gray_scale::Union{ClosedInterval,U} = 1000,
    calibration_value::Real = zero(T),
    background::Real = -1000,
    kwargs...,
) where {T<:Real,U<:Real}
    gray_scale = gray_scale isa Real ? gray_scale..gray_scale : gray_scale
    factor = 31 / 50 # magic number!
    #=
        If width or height are present, use them, otherwise use rows and cols.
        Finally, if rows or cols are not present, compute defaults.
    =#
    rows, cols = maybe(rows, height), maybe(cols, width)
    width = maybe(500, cols)
    height = isnothing(rows) ? round(Int, width * factor) : rows
    circ_zoom = 1
    swidth = round(Int, 17 / 25 * width)
    sheight = round(Int, 6 / 31 * height)
    pad = round(Int, 3 / 25 * width / circ_zoom)
    dist = 2pad ÷ 3
    radius = isnothing(radius) ? round(Int, width * circ_zoom / 17) : radius
    rows, cols = height, width
    imp = ImageParams(
        T;
        rows,
        cols,
        swidth,
        sheight,
        pad,
        dist,
        radius,
        gray_scale,
        calibration_value,
        background,
    )
    WhiteRect(imp; kwargs...)
end


function WhiteRect(geometry::AbstractParallelBeamGeometry; kwargs...)
    WhiteRect(
        datatype(geometry);
        geometry.width,
        geometry.height,
        kwargs...
    )
end


const _gs_images = (:GrayScaleLine, :GrayScalePyramid, :WhiteRect, :CircleImage)

for nm ∈ _gs_images
    @eval datatype(gs::$nm{T}) where T = T
end


for nm ∈ Geometry._geometry_names
    @eval begin
        $nm(imp::ImageParams; kwargs...) =
            $nm(datatype(imp); imp.width, imp.height, kwargs...)
        $nm(gs::AbstractTestImage; kwargs...) =
            $nm(datatype(gs); gs.width, gs.height, kwargs...)
    end
end


plateau_length(grsc::GrayScalePyramid) = round(Int, grsc.swidth * grsc.plateau)


project_image(image::AbstractTestImage, alg::AbstractProjectionAlgorithm; kwargs...) =
    project_image(image.image, alg; kwargs...)


for nm ∈ (:calibrate_image, :calibrate_tomogram)
    @eval begin
        $nm(imp::AbstractImageParams; kwargs...) =
            x -> $nm(x, imp; kwargs...)
        $nm(x, grsc::AbstractTestImage; kwargs...) =
            $nm(x, grsc.params; kwargs...)
        $nm(grsc::AbstractTestImage; kwargs...) =
            $nm(grsc.params; kwargs...)
        function $nm(
            m::Maybe,
            x::M;
            kwargs...
        ) where {M <: Union{AbstractImageParams,AbstractTestImage}}
            m ↣ $nm(x; kwargs...)
        end
    end
end


function calibration_data(imp::AbstractImageParams)
    min_pos = background_position(imp)
    max_pos = circle_position(imp)
    min_pos, max_pos
end


function rescale!(img::AbstractTestImage, args...; kwargs...)
    rescale!(img.image, args...; kwargs...)
end


function rescale(img::AbstractTestImage, args...; kwargs...)
    rescale(img.image, args...; kwargs...)
end

end # module
