module CTPlots

Base.Experimental.@optlevel 1

export plot_image, plot_sinogram, plot_tomogram
export plot_gray_scale

using RecipesBase: plot, plot!
using RecipesBase
using ..Monads
using ..CTImages: CTImage, CTSinogram, CTTomogram, ctfn, CTImageOrTomog
using ..CTImages: ctimage, ctsinogram, cttomogram
using ..TestImages: AbstractTestImage, ImageParams, gray_scale_indices
using ..CTScan: AbstractCTScanner
using ..Utils: linspace, ORI
using IntervalSets

function __init__()
    ENV["GKS_ENCODING"] = "utf-8"
end


"""
    plot_image(image::AbstractMatrix{T}; <keyword arguments>) where {T}

Plot CT image `image` with predefined options.

Additional keyword arguments are passed to the Plots library functions.

See also: [`plot_sinogram`](@ref), [`plot_tomogram`](@ref)
"""
function plot_image(image::AbstractMatrix{T}; kwargs...) where {T}
    plot(CTImage(image))
end


"""
    plot_tomogram(tomog::AbstractMatrix{T}; <keyword arguments>) where {T}

Plot reconstructed CT image `tomog` with predefined options.

Alias for `plot_image`. See [`plot_image`](@ref) for full reference.
Additional keyword arguments are passed to the Plots library.

See also: [`plot_sinogram`](@ref), [`plot_image`](@ref)
"""
function plot_tomogram(tomog::AbstractMatrix{T}; kwargs...) where {T}
    plot(CTTomogram(tomog); kwargs...)
end


"""
    plot_sinogram(sinog::AbstractMatrix{T}; <keyword arguments>) where {T}

Plot sinogram `sinog` with predefined options.

Additional keyword arguments are passed to the underlying plotting driver.

See also: [`plot_image`](@ref), [`plot_tomogram`](@ref)
"""
function plot_sinogram(sinog::AbstractMatrix{T}; kwargs...) where {T}
    plot(CTSinogram(sinog); kwargs...)
end

# TODO: Check lims for plot again!

@recipe function f(
    xs::AbstractVector, ys::AbstractVector, image::CTImageOrTomog
)
    rows, cols = size(image)
    proj = get(plotattributes, :projection, nothing)
    seriescolor --> :grays
    aspect_ratio --> :equal
    seriestype --> :heatmap
    _seriestype = get(plotattributes, :seriestype, nothing)
    if proj === :polar
        yaxis --> false
    elseif _seriestype === :heatmap
        tick_direction --> :out
        # xlims --> ((xs[1], xs[end]) .+ (-0.5, 0.5))
        xlims --> (xs[1], xs[end])
        # ylims --> ((ys[1], ys[end]) .+ (-0.5, 0.5))
        ylims --> (ys[1], ys[end])
    end
    xs, ys, mjoin(image)
end


@recipe function f(
    image::CTImageOrTomog,
    α::Optional{Real} = nothing,
    β::Optional{Real} = nothing,
)
    rows, cols = size(image)
    proj = get(plotattributes, :projection, nothing)
    if proj === :polar
        α = maybe(rows, α)
        β = maybe(2π, β)
        rs = linspace(0..α, rows+1)
        ϕs = linspace(0..β, cols+1)
        ϕs, rs, image
    else
        β = maybe(maybe((rows - 1) / 2, α), β)
        α = maybe((cols - 1) / 2, α)
        xs = linspace(-α..α, cols)
        ys = linspace(-β..β, rows)
        xs, ys, image
    end
end


const _sinog_xticks = [45i for i in 0:8]


@recipe function f(sinog::CTSinogram, α::Optional{Real} = nothing)
    nd, nϕ = size(sinog)
    xs = linspace(ORI(0..360), nϕ) # Now we use midpoints and not edges to support more backends
    α = maybe((nd - 1) / 2, α)
    ys = linspace(-α..α, nd)
    seriestype --> :heatmap
    seriescolor --> :grays
    xticks --> _sinog_xticks
    _seriestype = get(plotattributes, :seriestype, nothing)
    if _seriestype === :heatmap
        xrotation --> -45
        tick_direction --> :out
        # xlims --> ((xs[1], xs[end]) .+ (-0.5, 0.5) .* Δϕ)
        xlims --> (xs[1], xs[end])
        # ylims --> ((ys[1], ys[end]) .+ (-0.5, 0.5))
        ylims --> (ys[1], ys[end])
    end
    xs, ys, mjoin(sinog)
end


@recipe function f(
    xs::AbstractVector,
    ys::AbstractVector,
    gs::AbstractTestImage
)
    xs, ys, gs.image
end

@recipe function f(
    gs::AbstractTestImage,
    α::Optional{Real} = nothing,
    β::Optional{Real} = nothing
)
    gs.image, α, β
end


@recipe function f(imp::ImageParams, gray_scale_data::AbstractVector...)
    _, xs = gray_scale_indices(imp)
    seriestype --> :scatter
    linewidth --> 1.3
    markersize --> 1.8
    markerstrokewidth --> 0
    legend --> :outertopright
    x --> (xs .- imp.cols ÷ 2)
    gray_scale_data
end


@recipe function f(gst::AbstractCTScanner, s::Symbol)
    gst, Val(s)
end


for (nm, ctf) ∈ pairs(ctfn)
    v = Val{nm}
    @eval begin
        @recipe function f(gst::AbstractCTScanner, ::$v)
            mjoin($ctf(gst))
        end
    end
end


"""
    plot_image(gst::AbstractCTScanner)

Plot input image created for the test `gst`.
"""
plot_image(gst::AbstractCTScanner; kwargs...) = plot(gst, :image; kwargs...)


"""
    plot_sinogram(gst::AbstractCTScanner)

Plot the sinogram computed for the test `gst`.
"""
plot_sinogram(gst::AbstractCTScanner; kwargs...) = plot(gst, :sinog; kwargs...)


"""
    plot_tomogram(gst::AbstractCTScanner)

Plot the reconstructed image for the test `gst`.
"""
plot_tomogram(gst::AbstractCTScanner; kwargs...) = plot(gst, :tomog; kwargs...)


function _get_options(options, i)
    if isnothing(options)
        return (;)
    end
    if length(options) == 0
        return (;)
    end
    if i > length(options)
        return (;)
    end
    options[i]
end


"""
    plot_gray_scale(
        imp::ImageParams,
        gray_scale_data...;
        options = nothing,
        <keyword arguments>
    )

Plot gray scale data.

Additional keyword arguments are passed to the Plots library.
"""
function plot_gray_scale(
    imp::ImageParams,
    data::AbstractVector...;
    options = nothing,
    kwargs...)
    p = plot(imp, first(data); st = :line)
    for i ∈ 2:length(data)
        opts = _get_options(options, i)
        plot!(p, imp, data[i]; opts..., kwargs...)
    end
    p
end

end # module
