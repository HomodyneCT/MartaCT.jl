module TestImages

export indicator_matrix, indicator_matrix_2
export AbstractImageParams, ImageParams
export circle_position, background_position
export gray_scale_indices
export circle_image, gray_scale_image
export square_image, circle_polar_image
export combined_images_size, combine_images
export create_image, plateau_length
export AbstractGrayScale, GrayScaleLine, GrayScalePyramid
export WhiteRect


using IntervalSets
using ..Monads
using ..CTImages: CTImage
using ..Geometry: AbstractParallelBeamGeometry
using ..AbstractAlgorithms: AbstractProjectionAlgorithm
using ..Marta: datatype
import ..CTIO: yaml_repr, struct2dict
import ..Info: CTInfo
import ..AbstractAlgorithms: project_image
import ..CalibrationBase: calibrate_image, calibrate_tomogram
import Base: getproperty, size


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


abstract type AbstractImageParams end

"""
    struct ImageParams{T<:Real}

Hold information to construct the input test image.
"""
struct ImageParams{T<:Real} <: AbstractImageParams
    width::Int # Rectangle scale width
    height::Int # Rectangle scale height
    pad::Int # Padding around images
    dist::Int # Distance between images
    radius::Int # Calibration circle radius
    rows::Int # Test image rows
    cols::Int # Test image cols
    gray_scale::ClosedInterval{T} # Gray scale interval.
    calibration_value::T # Value for algorithm calibration.
    background::T

    function ImageParams{T}(
        width,
        height,
        pad,
        dist,
        radius,
        gray_scale::ClosedInterval = 0..1,
        calibration_value = nothing,
        background = nothing,
    ) where {T<:Real}
        rows, cols = combined_images_size(
            width = width,
            height = height,
            radius = radius,
            pad = pad,
            dist = dist,
        )

        minv, maxv = endpoints(gray_scale)

        if isnothing(background)
            background = minv
        end

        if isnothing(calibration_value)
            calibration_value = maxv
        end

        new(
            width,
            height,
            pad,
            dist,
            radius,
            rows,
            cols,
            gray_scale,
            T(calibration_value),
            T(background),
        )
    end
end


yaml_repr(imp::ImageParams) = struct2dict(imp)

function CTInfo(imp::ImageParams)
    CTInfo(pairs(struct2dict(imp))...)
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
    width = 200,
    height = 40,
    pad = 30,
    dist = 10,
    radius = 15,
    gray_scale::ClosedInterval = 0..1,
    calibration_value = nothing,
    background = nothing,
) where {T<:Real}
    ImageParams{T}(
        width,
        height,
        pad,
        dist,
        radius,
        gray_scale,
        calibration_value,
        background,
    )
end


show(io::IO, p::ImageParams) = print(
    io,
    """Image Params:
- width: $(p.width)
- height: $(p.height)
- pad: $(p.pad)
- dist: $(p.dist)
- radius: $(p.radius)
- size (W×H): $(p.cols) × $(p.rows)
- gray scale: $(p.gray_scale)
- calibration value: $(p.calibration_value)
- background value: $(p.background)""",
)


"""
    circle_position(imp::ImageParams)

Return circle position inside image given image parameters.
"""
circle_position(imp::ImageParams) = imp.pad + imp.radius + 1, imp.cols ÷ 2


"""
    background_position(imp::ImageParams)

Return suitable position to calibrate background.
"""
function background_position(imp::ImageParams)
    cr, cc = circle_position(imp)
    cr, cc ÷ 2
end


"""
    gray_scale_indices(imp::ImageParams)

Return a tuple `(row, column_range)` where `row` is the row index of the scale and `column_range` is the range of columns.
"""
function gray_scale_indices(imp::ImageParams)
    # row = imp.pad + 2 * imp.radius + imp.dist + imp.height ÷ 2 + 1
    min_row = imp.pad + 2 * imp.radius + imp.dist + 1
    max_row = min_row + imp.height - 1
    min_col = imp.pad + 1
    max_col = min_col + imp.width - 1
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
    radius = 20,
    value = 1,
    background = 0,
) where {T<:Real}
    r2 = radius^2
    rows = cols = 2radius
    circle_data = Matrix{T}(undef, rows, cols)
    for j in 1:cols, i in 1:rows
        x, y = j - radius, i - radius
        circle_data[i, j] = T(background)
        if x^2 + y^2 <= r2
            circle_data[i, j] = T(value)
        end
    end
    circle_data
end


"""
    circle_image(imp::ImageParams)

Create a square image with a circle of given value from parameters `imp`.
"""
function circle_image(imp::ImageParams{T}) where {T<:Real}
    circle_image(
        T;
        radius = imp.radius,
        value = imp.calibration_value,
        background = imp.background,
    )
end


"""
    gray_scale_image([T=Float32]; width=200, height=40, gray_scale=(0,1)) where {T <: Real}

Create an image with a gray scale rectangle with given scale gray_scale.

# Examples
```julia-repl
julia> gray_scale_image(width=80, height=40, gray_scale=(0,1))
40×80 Array{Float32,2}:
[...]
```

See also: [`circle_image`](@ref), [`combine_images`](@ref)
"""
function gray_scale_image(
    ::Type{T} = Float32;
    width = 200,
    height = 40,
    gray_scale::ClosedInterval = 0..1,
    background = nothing,
) where {T<:Real}
    minv, maxv = endpoints(gray_scale)

    if isnothing(background)
        background = minv
    end

    if minv == maxv
        val_range = fill(minv, width)
    else
        val_range = range(minv; stop = maxv, length = width) |> collect
    end

    rows, cols = height, width

    image = Matrix{T}(undef, rows, cols)
    image .= background

    start_r, end_r = 1, height
    start_c, end_c = 1, width

    for j in start_c:end_c, i in start_r:end_r
        image[i, j] = val_range[j]
    end

    image
end


"""
    pyramid_gray_scale_image([T=Float32]; width=200, height=40, gray_scale=(0,1)plateau = 0) where {T <: Real}

Create an image with a pyramid gray scale rectangle with given scale gray_scale.

# Examples
```julia-repl
julia> pyramid_gray_scale_image(width=80, height=40, gray_scale=(0,1))
40×80 Array{Float32,2}:
[...]
```

See also: [`gray_scale_image`](@ref), [`circle_image`](@ref), [`combine_images`](@ref)
"""
function pyramid_gray_scale_image(
    ::Type{T} = Float32;
    width = 200,
    height = 40,
    gray_scale::ClosedInterval = 0..1,
    background = nothing,
    plateau = 0,
) where {T<:Real}
    @assert 0 ≤ plateau ≤ 1

    minv, maxv = endpoints(gray_scale)

    if isnothing(background)
        background = minv
    end

    val_range = fill(T(maxv), width)

    plateau_len = round(Int, width * plateau)
    len = (width - plateau_len) ÷ 2

    if len ≠ 0
        half_scale = range(
            T(minv);
            step = T((maxv - minv) / (width - plateau_len - len - 1)),
            length = len,
        )

        val_range[1:len] .= half_scale
        val_range[end:-1:(width-len+1)] .= half_scale # Odd case included.
    end

    image = Matrix{T}(undef, height, width)
    image .= background

    start_r, end_r = 1, height
    start_c, end_c = 1, width

    for j in start_c:end_c, i in start_r:end_r
        image[i, j] = val_range[j]
    end

    image
end


"""
    gray_scale_image(imp::ImageParams)

Create an image with a gray scale rectangle with given scale gray_scale
from parameters `imp`.
"""
function gray_scale_image(imp::ImageParams{T}) where {T<:Real}
    gray_scale_image(T; imp.width, imp.height, imp.gray_scale, imp.background)
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
    plateau = 0,
) where {T<:Real}
    pyramid_gray_scale_image(
        T;
        imp.width,
        imp.height,
        imp.gray_scale,
        imp.background,
        plateau,
    )
end


function combined_images_size(; width, height, radius, pad, dist)
    height + 2radius + 2pad + dist, width + 2pad
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

    #@assert (rect_rows > circ_rows || rect_cols > circ_cols) "Circle image should be smaller than gray scale rectangle"

    radius = circ_rows ÷ 2

    rows, cols = combined_images_size(
        width = rect_cols,
        height = rect_rows,
        radius = radius,
        pad = imp.pad,
        dist = imp.dist,
    )

    image = Matrix{eltype(rect)}(undef, rows, cols)
    image .= imp.background

    crow_beg = imp.pad + 1
    crow_end = crow_beg + circ_rows - 1
    ccol_beg = cols ÷ 2 - radius + 1
    ccol_end = ccol_beg + circ_cols - 1

    image[crow_beg:crow_end, ccol_beg:ccol_end] .= circle

    rrow_beg = crow_end + imp.dist + 1
    rrow_end = rrow_beg + rect_rows - 1
    rcol_beg = imp.pad + 1
    rcol_end = rcol_beg + rect_cols - 1

    image[rrow_beg:rrow_end, rcol_beg:rcol_end] .= rect

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
function create_pyramid_image(par::ImageParams; plateau = 0)
    circ_img = circle_image(par)
    rect_img = pyramid_gray_scale_image(par; plateau)
    combine_images(par, rect_img, circ_img) |> CTImage
end


abstract type AbstractGrayScale end

size(gs::AbstractGrayScale) = size(gs.image)
size(gs::AbstractGrayScale, d::Integer) = size(gs.image, d)

struct GrayScaleLine{T<:Real} <: AbstractGrayScale
    params::ImageParams{T}
    image::CTImage{Matrix{T}}
end


function getproperty(grimg::GrayScaleLine, s::Symbol)
    s ∈ fieldnames(GrayScaleLine) && return getfield(grimg, s)
    s ≡ :width && return grimg.cols
    s ≡ :height && return grimg.rows
    getfield(grimg.params, s)
end


function GrayScaleLine(
    ::Type{T} = Float32;
    width = 500,
    height = nothing,
    gray_scale::ClosedInterval = -1000..1000,
    calibration_value = 0,
    background = -1000,
) where {T<:Real}
    factor = 31 / 50 # magic number!
    if isnothing(height)
        height = round(Int, width * factor)
    end
    gray_scale_width = round(Int, 17 / 25 * width)
    gray_scale_height = round(Int, 6 / 31 * height)
    pad = round(Int, 3 / 25 * width)
    dist = 2pad ÷ 3
    radius = round(Int, 1 / 17 * width)

    imp = ImageParams(
        T;
        width = gray_scale_width,
        height = gray_scale_height,
        pad,
        dist,
        radius,
        gray_scale,
        calibration_value,
        background,
    )

    GrayScaleLine(imp, create_image(imp))
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
    s ≡ :width && return grimg.cols
    s ≡ :height && return grimg.rows
    getfield(grimg.params, s)
end


function GrayScalePyramid(
    ::Type{T} = Float32;
    width = 500,
    height = nothing,
    gray_scale::ClosedInterval = -1000..1000,
    calibration_value = 0,
    background = -1000,
    plateau = 0,
) where {T<:Real}
    factor = 31 / 50 # magic number!
    if isnothing(height)
        height = round(Int, width * factor)
    end
    gray_scale_width = round(Int, 17 / 25 * width)
    gray_scale_height = round(Int, 6 / 31 * height)
    pad = round(Int, 3 / 25 * width)
    dist = 2pad ÷ 3
    radius = round(Int, 1 / 17 * width)

    imp = ImageParams(
        T;
        width = gray_scale_width,
        height = gray_scale_height,
        pad,
        dist,
        radius,
        gray_scale,
        calibration_value,
        background,
    )

    image = create_pyramid_image(imp; plateau)

    GrayScalePyramid(imp, T(plateau), image)
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
    nr,
    nϕ,
    radius;
    value = 1,
    background = 0,
) where {T<:Real}
    @assert radius <= nr "Radius should be less than number of radial points"
    mat = Matrix{T}(undef, nϕ, nr)
    mat .= background
    mat[:, 1:radius] .= T(value)
    mat
end
circle_polar_image(::Type{T}, nr, radius; value = 1) where {T} =
    circle_polar_image(T, nr, nr, radius; value = value)
circle_polar_image(nr, nϕ, radius; value = 1) =
    circle_polar_image(Float32, nr, nϕ, radius; value = value)
circle_polar_image(nr, radius; value = 1) =
    circle_polar_image(nr, nr, radius; value = value)
circle_polar_image(::Type{T}, radius; value = 1) where {T} =
    circle_polar_image(T, radius, radius; value = value)
circle_polar_image(radius; value = 1) where {T} =
    circle_polar_image(radius, radius; value = value)


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
function square_image(::Type{T}, r, c; l = nothing) where {T<:Real}
    matrix = zeros(T, r, c)

    if isnothing(l)
        l = min(r, c)
    end

    @assert l <= min(r, c)

    rin = r ÷ 2 - l ÷ 2 + 1
    rfi = rin + l - 1
    cin = c ÷ 2 - l ÷ 2 + 1
    cfi = cin + l - 1

    matrix[rin:rfi, cin:cfi] .= 1

    matrix
end


square_image(r, c; kwargs...) = square_image(Float32, r, c; kwargs...)


struct WhiteRect{T<:Real} <: AbstractGrayScale
    params::ImageParams{T}
    image::CTImage{Matrix{T}}
end


function getproperty(grimg::WhiteRect, s::Symbol)
    s ∈ fieldnames(WhiteRect) && return getfield(grimg, s)
    s ≡ :width && return grimg.cols
    s ≡ :height && return grimg.rows
    getfield(grimg.params, s)
end


function WhiteRect(
    ::Type{T} = Float32;
    width = 500,
    height = nothing,
    gray_scale::Union{ClosedInterval,U} = 1000,
    calibration_value = 0,
    background = -1000,
) where {T<:Real,U<:Real}
    if gray_scale isa Real
        gray_scale = gray_scale..gray_scale
    end

    grimg = GrayScaleLine(T; width, height, gray_scale, calibration_value, background)

    WhiteRect(grimg.params, grimg.image)
end


function WhiteRect(geometry::AbstractParallelBeamGeometry; kwargs...)
    WhiteRect(
        datatype(geometry);
        geometry.width,
        geometry.height,
        kwargs...
    )
end


plateau_length(grsc::GrayScalePyramid) = round(Int, grsc.width * grsc.plateau)


project_image(image::AbstractGrayScale, alg::AbstractProjectionAlgorithm; kwargs...) =
    project_image(image.image, alg; kwargs...)


for nm ∈ (:calibrate_image, :calibrate_tomogram)
    @eval begin
        $nm(imp::AbstractImageParams; kwargs...) =
            x -> $nm(x, imp; kwargs...)
        $nm(x, grsc::AbstractGrayScale; kwargs...) =
            $nm(x, grsc.params; kwargs...)
        $nm(grsc::AbstractGrayScale; kwargs...) =
            $nm(grsc.params; kwargs...)
        $nm(
            m::Maybe,
            x::Union{AbstractImageParams,AbstractGrayScale};
            kwargs...
        ) = m ↣ $nm(x; kwargs...)
    end
end

end # module
