module CTScan

Base.Experimental.@optlevel 0

export AbstractCTScanner, CTScanner, FBPScanner
export rename, transform_data
export compute_gray_scale
export make_phantom, create_image!
export project_image!, reconstruct_image!
export para2fan!, fan2para!
export compute_rmse, compute_χ2, fit_gray_scale
export load_image!, load_sinogram!, load_tomogram!
export calibrate_image!, calibrate_tomogram!
export sample_sinogram!

import Base: show, copy, similar, getproperty
using IntervalSets
import Statistics; const Stat = Statistics
import CurveFit; const Fit = CurveFit
using ..Applicative
using ..Monads
using ..Geometry
import ..AbstractAlgorithms: reconstruct_image, project_image
using ..AbstractAlgorithms
import ..Marta: datatype, _atype
using ..RadonAlgorithm
import ..CTImages: ctimage, ctsinogram, cttomogram
using ..CTImages
using ..CTData
import ..FanBeam: para2fan, fan2para
import ..CalibrationBase: calibrate_image, calibrate_tomogram
using ..Calibration
using ..TestImages
import ..CTIO: load_image, write_image
import ..CTIO: load_sinogram, write_sinogram
import ..CTIO: load_tomogram, write_tomogram


include("CTImageData.jl")
using .CTImageData: CTData


"""
    abstract type AbstractCTScanner{<:AbstractArray,<:AbstractGeometry} end

Abstract type for tests in this suite.
"""
abstract type AbstractCTScanner end

const CTScannerNameType = Union{Symbol,String}

function rename(gst::AbstractCTScanner, name::AbstractString)
    typeof(gst)(gst.geometry, gst.data; name)
end

_atype(gst::AbstractCTScanner) = _atype(gst.data)

ctimage(gst::AbstractCTScanner) = ctimage(gst.data)
ctsinogram(gst::AbstractCTScanner) = ctsinogram(gst.data)
cttomogram(gst::AbstractCTScanner) = cttomogram(gst.data)


similar(gst::AbstractCTScanner, img::AbstractCTImage...) = similar(gst, CTData(img...))

# function similar(
#     gst::AbstractCTScanner,
#     g::AbstractGeometry,
#     img::AbstractCTImage...,
# )
#     similar(gst, g, CTData(img...))
# end
function similar(gst::AbstractCTScanner, img::AbstractCTImage...)
    similar(gst, CTData(img...))
end


for (nm, ctf) ∈ pairs(CTImages.ctfn)
    v = Val{nm}
    @eval begin
        t = CTImages.ctnames.$nm
        function transform_data(f::Function, gst::AbstractCTScanner, ::$v)
            mbind($ctf(gst)) do img
                newdata = similar(gst.data, map(f, something(img)))
                similar(gst, newdata)
            end
        end
        transform_data(f::Function, gst::AbstractCTScanner, ::Type{t}) =
            transform_data(f, gst, $v())
    end
end


transform_data(f::Function, gst::AbstractCTScanner, s::Symbol) =
    transform_data(f, gst, Val(s))

struct CTScanner{
    Alg <: AbstractReconstructionAlgorithm,
    M <: AbstractArray{<:Real},
} <: AbstractCTScanner
    name::CTScannerNameType
    study_id::String
    algorithm::Alg
    data::CTData{M}

    """
    CTScanner{A,M}(
        alg::A[, data::CTData];
        name::Optional{Union{Symbol,String}} = A,
        study_id::Optional{String} = "Unknown"
    ) where {A<:AbstractReconstructionAlgorithm,M<:AbstractArray}

    Construct `CTScanner` from an algorithm.
    """
    function CTScanner{A,M}(
        alg::A,
        data::CTData{M} = CTData(M);
        name::Optional{CTScannerNameType} = nothing,
        study_id::Optional{String} = nothing,
    ) where {A <: AbstractReconstructionAlgorithm,M <: AbstractArray{<:Real}}
        new(maybe(nameof(A), name), maybe("Unknown", study_id), alg, data)
    end
end


const FBPScanner{M} = CTScanner{FBP,M}


# @inline function getproperty(ct::CTScanner, s::Symbol)
#     s ≡ :geometry && return ct.algorithm.geometry
#     if s ≡ :params
#         ct isa FBPScanner && return nothing
#     end
#     getfield(ct, s)
# end


function CTScanner(alg::AbstractReconstructionAlgorithm, data::CTData; kwargs...)
    CTScanner{typeof(alg), _atype(data)}(alg, data; kwargs...)
end


function CTScanner(alg::AbstractReconstructionAlgorithm, gst::AbstractCTScanner; kwargs...)
    CTScanner(alg, gst.data; kwargs...)
end


function similar(gst::CTScanner, data::Optional{<:CTData} = nothing)
    CTScanner(
        gst.algorithm,
        maybe(CTData(_atype(gst)), data);
        gst.name,
        gst.study_id,
    )
end


# """
#     FBPScanner(
#         geometry::AbstractGeometry,
#         data::CTData;
#         name::Optional{String},
#         study_id::Optional{String}
#     )

# Construct a `FBPScanner` object.
# """
# function FBPScanner(
#     geometry::G,
#     data::CTData{M} = CTData(datatype(G));
#     name::Optional{CTScannerNameType} = nothing,
#     study_id::Optional{String} = nothing,
# ) where {G<:AbstractGeometry,M<:AbstractArray}
#     FBPScanner{M}(FBP(geometry), data; name, study_id)
# end

"""
    FBPScanner(data::CTData; name::Optional{String}, study_id::Optional{String})

Construct a `FBPScanner` object.
"""
function FBPScanner(
    ::Type{T} = Float32,
    data::CTData{M} = CTData(T);
    name::Optional{CTScannerNameType} = nothing,
    study_id::Optional{String} = nothing,
) where {G<:AbstractGeometry,M<:AbstractArray}
    FBPScanner{M}(FBP(), data; name, study_id)
end


# """
#     FBPScanner(
#         geometry::AbstractGeometry,
#         image::AbstractGrayScale;
#         name::Optional{String},
#         study_id::Optional{String}
#     )

# Construct a `FBPScanner` object.
# """
# function FBPScanner(
#     geometry::AbstractGeometry,
#     gs::AbstractTestImage;
#     name::Optional{CTScannerNameType} = nothing,
#     study_id::Optional{String} = nothing,
# )
#     FBPScanner(geometry, gs.image; name, study_id)
# end
"""
    FBPScanner(
        image::AbstractGrayScale;
        name::Optional{String},
        study_id::Optional{String}
    )

Construct a `FBPScanner` object.
"""
function FBPScanner(
    gs::AbstractTestImage;
    name::Optional{CTScannerNameType} = nothing,
    study_id::Optional{String} = nothing,
)
    FBPScanner(gs.image; name, study_id)
end


# """
#     FBPScanner(
#         image::AbstractGrayScale;
#         name::Optional{String},
#         study_id::Optional{String}
#     )

# Construct a `FBPScanner` object.
# """
# function FBPScanner(
#     gs::AbstractTestImage;
#     name::Optional{CTScannerNameType} = nothing,
#     study_id::Optional{String} = nothing,
# )
#     FBPScanner(ParallelBeamGeometry(gs), gs; name, study_id)
# end


"""
    FBPScanner(gst::AbstractCTScanner; name::Optional{String}, study_id::Optional{String})

Construct a `FBPScanner` object from another `AbstractCTScanner` object. This allows to reuse data.
"""
function FBPScanner(
    gst::AbstractCTScanner;
    name::Optional{CTScannerNameType} = nothing,
    study_id::Optional{String} = nothing,
)
    data = CTData(ctimage(gst), ctsinogram(gst), nothing)
    name = maybe(gst.name, name)
    study_id = maybe(gst.study_id, study_id)
    #FBPScanner(gst.geometry, data; name, study_id)
    FBPScanner(eltype(data), data; name, study_id)
end


# function FBPScanner(g::AbstractGeometry, img::AbstractCTImage...; kwargs...)
#     FBPScanner(g, CTData(img...); kwargs...)
# end
function FBPScanner(img::AbstractCTImage...; kwargs...)
    FBPScanner(CTData(img...); kwargs...)
end


# function copy(gst::FBPScanner)
#     FBPScanner(gst.geometry, copy(gst.data); gst.name, gst.study_id)
# end
function copy(gst::FBPScanner)
    FBPScanner(copy(gst.data); gst.name, gst.study_id)
end


# function similar(
#     gst::FBPScanner,
#     g::AbstractGeometry,
#     data::Optional{CTData} = nothing,
# )
#     FBPScanner(
#         g,
#         maybe(gst.data, data);
#         gst.name,
#         gst.study_id,
#     )
# end
function similar(gst::FBPScanner, data::Optional{CTData} = nothing)
    FBPScanner(maybe(gst.data, data); gst.name, gst.study_id)
end


function show(io::IO, gst::CTScanner)
    print(
        io,
        """*** $(gst.name) Scanner ***
        $(gst.algorithm)""",
    )
end


"""
    create_image!(gst::AbstractCTScanner, par::ImageParams)

Create gray scale image for the test `gst`.
"""
function create_image!(gst::AbstractCTScanner, par::ImageParams)
    #similar(gst, similar(gst.data, create_image(par)))
    gst.data.image = create_image(par)
    gst
end


# """
#     project_image(gst::AbstractCTScanner, alg::Optional{A}; <keyword arguments>)
#         where {A <: AbstractProjectionAlgorithm}

# Compute sinogram for the test `gst`.

# Keyword arguments depend on the algorithm employed, please see the relative
# documentation.
# """
# function project_image(
#     gst::AbstractCTScanner,
#     alg::Optional{A} = nothing;
#     kwargs...
# ) where {A <: AbstractProjectionAlgorithm}
#     sinog = project_image(
#         ctimage(gst),
#         isnothing(alg) ? Radon(gst.geometry) : alg;
#         kwargs...
#     )
#     similar(gst, similar(gst.data, sinog))
# end

"""
project_image(gst::AbstractCTScanner, g::AbstractGeometry,
alg::Optional{A}; <keyword arguments>)
        where {A <: AbstractProjectionAlgorithm}

Compute sinogram for the test `gst`.

Keyword arguments depend on the algorithm employed, please see the relative
documentation.
"""
function project_image(
    gst::AbstractCTScanner,
    g::AbstractGeometry,
    alg::Optional{A} = nothing;
    kwargs...
) where {A <: AbstractProjectionAlgorithm}
    sinog = project_image(
        ctimage(gst),
        g,
        isnothing(alg) ? Radon() : alg;
        kwargs...
    )
    similar(gst, similar(gst.data, sinog))
end

function project_image!(
    gst::AbstractCTScanner,
    g::AbstractGeometry,
    alg::Optional{A} = nothing;
    kwargs...
) where {A <: AbstractProjectionAlgorithm}
    sinog = project_image(
        ctimage(gst),
        g,
        isnothing(alg) ? Radon() : alg;
        kwargs...
    )
    gst.data.sinog = sinog
    gst
end


# function para2fan(gst::AbstractCTScanner; kwargs...)
#     fbg, sinog_fan = ctsinogram(gst) ↣ para2fan(gst.geometry; kwargs...)
#     similar(gst, fbg, similar(gst.data, sinog_fan))
# end
function para2fan(gst::AbstractCTScanner, g::AbstractGeometry; kwargs...)
    fbg, sinog_fan = ctsinogram(gst) ↣ para2fan(g; kwargs...)
    similar(gst, fbg, similar(gst.data, sinog_fan))
end

function para2fan!(gst::AbstractCTScanner, g::AbstractGeometry; kwargs...)
    fbg, sinog_fan = ctsinogram(gst) ↣ para2fan(g; kwargs...)
    gst.data.sinog = sinog_fan
    gst
end


# function para2fan(gst::AbstractCTScanner, fbg::FanBeamGeometry; kwargs...)
#     fbg, sinog_fan = ctsinogram(gst) ↣ para2fan(fbg; kwargs...)
#     similar(gst, fbg, similar(gst.data, sinog_fan))
# end


# function fan2para(gst::AbstractCTScanner; kwargs...)
#     pbg, sinog_para = ctsinogram(gst) ↣ fan2para(gst.geometry; kwargs...)
#     similar(gst, pbg, similar(gst.data, sinog_para))
# end
function fan2para(gst::AbstractCTScanner; kwargs...)
    pbg, sinog_para = ctsinogram(gst) ↣ fan2para(gst.geometry; kwargs...)
    similar(gst, pbg, similar(gst.data, sinog_para))
end

function fan2para!(gst::AbstractCTScanner, g::AbstractGeometry; kwargs...)
    pbg, sinog_para = ctsinogram(gst) ↣ fan2para(g; kwargs...)
    gst.data.sinog = sinog_para
    gst
end


# """
#     reconstruct_image(gst::AbstractCTScanner; <keyword arguments>)

# Reconstruct CT image from test `gst`.

# Keyword arguments depend on the specific algorithm used, please see the relative
# documentation.
# """
# function reconstruct_image(gst::AbstractCTScanner; kwargs...)
#     tomog = reconstruct_image(ctsinogram(gst), gst.algorithm; kwargs...)
#     similar(gst, similar(gst.data, tomog))
# end

"""
    reconstruct_image(gst::AbstractCTScanner, g::AbstractGeometry; <keyword arguments>)

Reconstruct CT image from test `gst`.

Keyword arguments depend on the specific algorithm used, please see the relative
documentation.
"""
function reconstruct_image(gst::AbstractCTScanner, g::AbstractGeometry; kwargs...)
    tomog = reconstruct_image(ctsinogram(gst), g, gst.algorithm; kwargs...)
    similar(gst, similar(gst.data, tomog))
end

function reconstruct_image!(gst::AbstractCTScanner, g::AbstractGeometry; kwargs...)
    tomog = reconstruct_image(ctsinogram(gst), g, gst.algorithm; kwargs...)
    gst.data.tomog = tomog
    gst
end


"""
    load_sinogram(io::IO, gst::AbstractCTScanner)

Load sinogram from stream `io`.
"""
function load_sinogram(io::IO, gst::AbstractCTScanner)
    sinog = load_sinogram(io; gst.geometry.ns, gst.geometry.nϕ)
    similar(gst, similar(gst.data, sinog))
end

function load_sinogram!(io::IO, gst::AbstractCTScanner)
    sinog = load_sinogram(io; gst.geometry.ns, gst.geometry.nϕ)
    gst.data.sinog = sinog
    gst
end


"""
    load_sinogram(f::AbstractString, gst::AbstractCTScanner)

Load sinogram from file `f`.
"""
function load_sinogram(f::AbstractString, gst::AbstractCTScanner)
    open(f) do s
        load_sinogram(s, gst)
    end
end

function load_sinogram!(f::AbstractString, gst::AbstractCTScanner)
    open(f) do s
        load_sinogram!(s, gst)
    end
end


"""
    load_image(io::IO, gst::AbstractCTScanner)

Load image from stream `io`.
"""
function load_image(io::IO, gst::AbstractCTScanner)
    image = load_image(io; gst.geometry.rows, gst.geometry.cols)
    similar(gst, similar(gst.data, image))
end

function load_image!(io::IO, gst::AbstractCTScanner)
    image = load_image(io; gst.geometry.rows, gst.geometry.cols)
    gst.data.image = image
    gst
end


"""
    load_image(f::AbstractString, gst::AbstractCTScanner)

Load image from file `f`.
"""
function load_image(f::AbstractString, gst::AbstractCTScanner)
    open(f) do s
        load_image(s, gst)
    end
end

function load_image!(f::AbstractString, gst::AbstractCTScanner)
    open(f) do s
        load_image!(s, gst)
    end
end


"""
    write_image(io::IO, gst::AbstractCTScanner)

Write input image to stream `io`.
"""
function write_image(io::IO, gst::AbstractCTScanner)
    mbind(ctimage(gst)) do x
        write_image(io, x; header = false, row_major = true)
    end
end


"""
    write_image(f::AbstractString, gst::AbstractCTScanner)

Write input image to file `f`.
"""
function write_image(f::AbstractString, gst::AbstractCTScanner)
    open(f, "w") do s
        write_image(s, gst)
    end
end


"""
    write_sinogram(io::IO, gst::AbstractCTScanner)

Write sinogram to stream `io`.
"""
function write_sinogram(io::IO, gst::AbstractCTScanner)
    mbind(ctsinogram(gst)) do x
        write_sinogram(io, x; header = false)
    end
end


"""
    write_sinogram(f::AbstractString, gst::AbstractCTScanner)

Write sinogram to file `f`.
"""
function write_sinogram(f::AbstractString, gst::AbstractCTScanner)
    open(f, "w") do s
        write_sinogram(s, gst)
    end
end


"""
    load_tomogram(io::IO, gst::AbstractCTScanner)

Load tomogram from stream `io`.
"""
function load_tomogram(io::IO, gst::AbstractCTScanner)
    tomog = load_tomogram(io; gst.geometry.rows, gst.geometry.cols)
    similar(gst, similar(gst.data, tomog))
end

function load_tomogram!(io::IO, gst::AbstractCTScanner)
    tomog = load_tomogram(io; gst.geometry.rows, gst.geometry.cols)
    gst.data.tomog = tomog
    gst
end


"""
    load_tomogram(f::AbstractString, gst::AbstractCTScanner)

Load tomogram from file `f`.
"""
function load_tomogram(f::AbstractString, gst::AbstractCTScanner)
    open(f) do s
        load_tomogram(s, gst)
    end
end

function load_tomogram!(f::AbstractString, gst::AbstractCTScanner)
    open(f) do s
        load_tomogram!(s, gst)
    end
end


"""
    write_tomogram(io::IO, gst::AbstractCTScanner)

Write reconstructed image to stream `io`.
"""
function write_tomogram(io::IO, gst::AbstractCTScanner)
    mbind(cttomogram(gst)) do x
        write_tomogram(io, x; header = false, row_major = true)
    end
end


"""
    write_tomogram(f::AbstractString, gst::AbstractCTScanner)

Write tomogram to file `f`.
"""
function write_tomogram(f::AbstractString, gst::AbstractCTScanner)
    open(f, "w") do s
        write_tomogram(s, gst)
    end
end


"""
    sample_sinogram(gst::AbstractCTScanner; <keyword arguments>)

Returns a new test object with the resampled sinogram.
"""
function sample_sinogram(gst::AbstractCTScanner; kwargs...)
    sinog = ctsinogram(gst) ↣ sample_sinogram(; kwargs...)
    similar(gst, similar(gst.data, sinog))
end

function sample_sinogram!(gst::AbstractCTScanner; kwargs...)
    sinog = ctsinogram(gst) ↣ sample_sinogram(; kwargs...)
    gst.data.sinog = sinog
    gst
end


"""
    calibrate_image(gst::AbstractCTScanner, imp::ImageParams; interval=nothing, window=nothing)

Perform calibration of input image using image parameters as reference.
"""
function calibrate_image(
    gst::AbstractCTScanner,
    imp::AbstractImageParams;
    interval::Optional{ClosedInterval{T}} = nothing,
    window::Optional{ClosedInterval{U}} = nothing,
) where {T<:Real,U<:Real}
    cimg = calibrate_image(
        ctimage(gst),
        imp;
        interval,
        window,
    )
    similar(gst, similar(gst.data, cimg))
end

function calibrate_image!(
    gst::AbstractCTScanner,
    imp::AbstractImageParams;
    interval::Optional{ClosedInterval{T}} = nothing,
    window::Optional{ClosedInterval{U}} = nothing,
) where {T<:Real,U<:Real}
    cimg = calibrate_image(
        ctimage(gst),
        imp;
        interval,
        window,
    )
    gst.data.image = cimg
    gst
end


"""
    calibrate_image(gst::AbstractCTScanner; min_pos, max_pos, interval=0..1, window=nothing)

Perform calibration of input image with reference values.
"""
function calibrate_image(
    gst::AbstractCTScanner;
    min_pos,
    max_pos,
    interval::ClosedInterval = 0..1,
    window::Optional{ClosedInterval{T}} = nothing,
) where {T<:Real}
    cimg = calibrate_image(
        ctimage(gst);
        interval,
        min_pos,
        max_pos,
        window,
    )
    similar(gst, similar(gst.data, cimg))
end

function calibrate_image!(
    gst::AbstractCTScanner;
    min_pos,
    max_pos,
    interval::ClosedInterval = 0..1,
    window::Optional{ClosedInterval{T}} = nothing,
) where {T<:Real}
    cimg = calibrate_image(
        ctimage(gst);
        interval,
        min_pos,
        max_pos,
        window,
    )
    gst.data.image = cimg
    gst
end


"""
    calibrate_tomogram(gst::AbstractCTScanner, imp::ImageParams; interval=nothing, window=nothing)

Perform calibration of reconstructed image using image parameters as reference.
"""
function calibrate_tomogram(
    gst::AbstractCTScanner,
    imp::AbstractImageParams;
    interval::Optional{ClosedInterval{T}} = nothing,
    window::Optional{ClosedInterval{U}} = nothing,
) where {T<:Real,U<:Real}
    cimg = calibrate_tomogram(
        cttomogram(gst),
        imp;
        interval,
        window,
    )
    similar(gst, similar(gst.data, cimg))
end

function calibrate_tomogram!(
    gst::AbstractCTScanner,
    imp::AbstractImageParams;
    interval::Optional{ClosedInterval{T}} = nothing,
    window::Optional{ClosedInterval{U}} = nothing,
) where {T<:Real,U<:Real}
    cimg = calibrate_tomogram(
        cttomogram(gst),
        imp;
        interval,
        window,
    )
    gst.data.tomog = cimg
    gst
end


"""
    calibrate_tomogram(gst::AbstractCTScanner; min_pos, max_pos, interval=0..1, window=nothing)

Perform calibration of reconstructed image with reference values.
"""
function calibrate_tomogram(
    gst::AbstractCTScanner;
    min_pos,
    max_pos,
    interval::ClosedInterval = 0..1,
    window::Optional{ClosedInterval{T}} = nothing,
) where {T<:Real}
    cimg = calibrate_tomogram(
        cttomogram(gst);
        min_pos,
        max_pos,
        interval,
        window,
    )
    similar(gst, similar(gst.data, cimg))
end

function calibrate_tomogram!(
    gst::AbstractCTScanner;
    min_pos,
    max_pos,
    interval::ClosedInterval = 0..1,
    window::Optional{ClosedInterval{T}} = nothing,
) where {T<:Real}
    cimg = calibrate_tomogram(
        cttomogram(gst);
        min_pos,
        max_pos,
        interval,
        window,
    )
    gst.data.tomog = cimg
    gst
end


"""
    compute_gray_scale(image::AbstractMatrix, imp::ImageParams; mean = false)

Compute the corresponding gray scale for the test image `image`.
"""
function compute_gray_scale(image::AbstractMatrix, imp::ImageParams; mean = false)
    mean && return Stat.mean(image[gray_scale_indices(imp)...], dims = 1) |> vec
    rind, cind = gray_scale_indices(imp)
    i = sum(extrema(rind)) ÷ 2
    image[i,cind]
end


# function compute_gray_scale(
#     image::CTImageOrTomog, params::ImageParams; kwargs...)
#     mbind(image) do x
#         compute_gray_scale(x, params; kwargs...)
#     end
# end

compute_gray_scale(grsc::AbstractGrayScale; kwargs...) =
    compute_gray_scale(grsc.image, grsc.params)

compute_gray_scale(imp::ImageParams; kwargs...) = x -> compute_gray_scale(x, imp; kwargs...)


"""
    compute_gray_scale(gst::AbstractCTScanner, imp::ImageParams; mean = false)

Compute the corresponding gray scale for the reconstructed image.
"""
function compute_gray_scale(gst::AbstractCTScanner, grsc::AbstractGrayScale; kwargs...)
    cttomogram(gst) ↣ compute_gray_scale(grsc.params; kwargs...)
end


function make_phantom(f::Function, g::AbstractParallelBeamGeometry)
    g |> f |> CTImage
end


function compute_rmse(
    grimg::AbstractGrayScale,
    expected::AbstractVector,
    measured::AbstractVector;
    relative = false,
    scale::Optional{ClosedInterval} = nothing,
)
    scale = maybe(grimg.gray_scale, scale)
    rmse = Stat.mean(abs2, measured - expected) |> sqrt
    relative && return rmse / width(scale)
    rmse
end


function _fit_line(ys)
    xs = eachindex(ys)
    a, b = Fit.linear_fit(xs, ys)
    @. a + b * xs
end


function fit_gray_scale(::Union{GrayScaleLine,WhiteRect}, data::AbstractVector)
    _fit_line(data)
end


function fit_gray_scale(grimg::GrayScalePyramid, data::AbstractVector)
    plateau_len = plateau_length(grimg)
    len = (grimg.width - plateau_len) ÷ 2
    len == 0 && return _fit_line(data)
    steps = accumulate((len, len + plateau_len, grimg.width), init = 0) do x, y
        (last(x) + 1):y
    end
    lines = _fit_line.(data[l] for l ∈ steps)
    reduce(vcat, lines)
end


function compute_rmse(grimg::Union{GrayScaleLine,WhiteRect}, data::AbstractVector; relative = false, scale = nothing)
    fitted_data = fit_gray_scale(grimg, data)
    compute_rmse(grimg, fitted_data, data; relative, scale)
end


function compute_rmse(grimg::GrayScalePyramid, data::AbstractVector; relative = false, scale = nothing)
    fitted_data = fit_gray_scale(grimg, data)
    compute_rmse(grimg, fitted_data, data; relative, scale)
end


function compute_χ2(grimg::AbstractGrayScale, expected::AbstractVector, measured::Vector, variances::Vector)
    ν = length(expected) - 2
    sum(@. (measured - expected) ^ 2 / variances) / ν
end

end # module
