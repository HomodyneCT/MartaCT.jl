function iradon_fast_threaded(
    sinog::AbstractMatrix{T};
    rows::Optional{Integer} = nothing,
    cols::Optional{Integer} = nothing,
    filter::Optional{F} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    interpolation::Optional{Interp} = nothing,
) where {
    T <: Real,
    F <: CTFilterOrFunc,
    Interp <: Union{Function,AbstractBilinearInterpolation}
}
    M = typeof(sinog)
    nd, nϕ = size(sinog)
    rows = maybe(round(Int, nd / sqrt(2)), rows)
    cols = maybe(rows, cols)

    interpolation = maybe(interpolate, interpolation)

    filtered = maybe(RamLak(), filter) do f
        filter_freq = fft(sinog, 1) .* f(T, nd, nϕ)
        ifft(filter_freq, 1) .|> real .|> T # fastest solution (?),
                                            # not '.|> T ∘ real'
    end

    interp = interpolation(filtered)

    Δϕ::T = deg2rad(α) / nϕ
    ϕ₀::T = deg2rad(α₀)
    x₀::T = T(cols - 1) / 2
    y₀::T = T(rows - 1) / 2
    t₀::T = T(nd - 1) / 2

    # ϕs = collect(enumerate((0:nϕ-1) * Δϕ .+ ϕ₀))
    ϕs = collect(enumerate(range(ϕ₀; step = Δϕ, length = nϕ)))
    xs = range(-1, 1; length = cols)
    # xs = range(-one(T); step = T(2 / cols), length = cols)
    ys = range(-1, 1; length = rows)
    # ys = range(-one(T); step = T(2/rows), length = rows)
    indices = Vector{Tuple{T,T}}(undef, rows * cols)

    for k in LinearIndices(indices)
        @inbounds indices[k] = xs[(k - 1) ÷ rows + 1], ys[(k - 1) % rows + 1]
    end

    @info "Computing inverse Radon transform..."
    p = Progress(nϕ, 0.2)

    temp_images = Vector{Matrix{T}}(undef, Threads.nthreads())
    fill!(temp_images, zeros(T, rows, cols))

    Threads.@threads for (iϕ, ϕ) in ϕs
        cϕ, sϕ = cos(ϕ), sin(ϕ)
        id = Threads.threadid()
        @inbounds img = temp_images[id]
        for (k, (x, y)) in enumerate(indices)
            # To be consistent with our conventions should be '+'.
            t = (x * cϕ + y * sϕ + 1) * t₀ + 1
            if 1 ≤ t ≤ nd
                @inbounds img[k] += interp(t, T(iϕ))
            end
        end
        next!(p)
    end

    temp_images |> sum |> CTTomogram{M}
end


function iradon_fast_threaded(
    sinog::AbstractMatrix{T},
    geometry::AbstractParallelBeamGeometry,
    filter::Optional{F} = nothing;
    kwargs...
) where {T <: Real, F <: CTFilterOrFunc}
    rows = geometry.rows
    cols = geometry.cols
    α = geometry.α
    α₀ = geometry.α₀
    iradon_fast_threaded(sinog; rows, cols, α, α₀, filter, kwargs...)
end


function iradon_fast_threaded(
    sinog::AbstractMatrix,
    filter::CTFilterOrFunc;
    geometry::Optional{G} = nothing,
    kwargs...
) where {G <: AbstractParallelBeamGeometry}
    isnothing(geometry) &&  return iradon_fast_threaded(sinog; filter, kwargs...)
    iradon_fast_threaded(sinog, geometry, filter; kwargs...)
end


iradon_fast_threaded(; kwargs...) = x -> iradon_fast_threaded(x; kwargs...)
iradon_fast_threaded(g::AbstractParallelBeamGeometry, f::CTFilterOrFunc; kwargs...) =
    x -> iradon_fast_threaded(x, g, f; kwargs...)
iradon_fast_threaded(g::AbstractParallelBeamGeometry; kwargs...) =
    x -> iradon_fast_threaded(x, g; kwargs...)
    iradon_fast_threaded(f::CTFilterOrFunc; kwargs...) =
    x -> iradon_fast_threaded(x, f; kwargs...)


function iradon_threaded(
    sinog::AbstractMatrix{T};
    rows::Optional{Integer} = nothing,
    cols::Optional{Integer} = nothing,
    filter::Optional{F} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    interpolation::Optional{Interp} = nothing,
) where {
    T <: Real,
    F <: CTFilterOrFunc,
    Interp <: Union{Function,AbstractBilinearInterpolation}
}
    M = typeof(sinog)
    nd, nϕ = size(sinog)
    rows = maybe(round(Int, nd / sqrt(2)), rows)
    cols = maybe(rows, cols)

    filtered = maybe(RamLak(), filter) do f
        filter_freq = fft(sinog, 1) .* f(T, nd, nϕ)
        ifft(filter_freq, 1) .|> real .|> T # fastest solution (?),
                                            # not '.|> T ∘ real'
    end

    interp = isnothing(interpolation) ? interpolate(filtered) : interpolation(filtered)

    Δϕ::T = deg2rad(α) / nϕ
    ϕ₀::T = deg2rad(α₀)
    x₀::T = T(cols - 1) / 2
    y₀::T = T(rows - 1) / 2
    t₀::T = T(nd - 1) / 2

    image = zeros(T, rows, cols)
    xs = range(-1, 1; length = cols)
    # xs = range(-1; step = 2 / cols, length = cols)
    ys = range(-1, 1; length = rows)
    # ys = range(-1; step = 2 / rows, length = rows)
    indices = Vector{Tuple{Int,T,T}}(undef, length(image))

    for k in LinearIndices(image)
        @inbounds indices[k] = k, xs[(k - 1) ÷ rows + 1], ys[(k - 1) % rows + 1]
    end

    # ϕs = map(enumerate((0:nϕ-1) * Δϕ .+ ϕ₀)) do (iϕ, ϕ)
    #     T(iϕ), cos(ϕ), sin(ϕ)
    # end
    indϕs = enumerate(range(ϕ₀; step = Δϕ, length = nϕ))
    csϕs = map(indϕs) do (iϕ, ϕ)
        T(iϕ), cos(ϕ), sin(ϕ)
    end

    p = Progress(length(image), 0.2, "Computing inverse Radon transform: ")

    Threads.@threads for (k, x, y) in indices
        @inbounds image[k] = sum(csϕs) do (iϕ, cϕ, sϕ)
            # To be consistent with our conventions should be '+'.
            t::T = (x * cϕ + y * sϕ + 1) * t₀ + 1
            1 ≤ t ≤ nd && return interp(t, iϕ)
        end
        next!(p)
    end

    image |> CTTomogram{M}
end


function iradon_threaded(
    sinog::AbstractMatrix{T},
    geometry::AbstractParallelBeamGeometry,
    filter::Optional{F} = nothing;
    kwargs...
) where {T <: Real, F <: CTFilterOrFunc}
    rows = geometry.rows
    cols = geometry.cols
    α = geometry.α
    α₀ = geometry.α₀
    iradon_threaded(sinog; rows, cols, α, α₀, filter, kwargs...)
end


function iradon_threaded(
    sinog::AbstractMatrix,
    filter::CTFilterOrFunc;
    geometry::Optional{G} = nothing,
    kwargs...
) where {G <: AbstractParallelBeamGeometry}
    isnothing(geometry) && return iradon_threaded(sinog; filter, kwargs...)
    iradon_threaded(sinog, geometry, filter; kwargs...)
end


iradon_threaded(; kwargs...) = x -> iradon_threaded(x; kwargs...)
iradon_threaded(g::AbstractParallelBeamGeometry, f::CTFilterOrFunc; kwargs...) =
    x -> iradon_threaded(x, g, f; kwargs...)
iradon_threaded(g::AbstractParallelBeamGeometry; kwargs...) =
    x -> iradon_threaded(x, g; kwargs...)
iradon_threaded(f::CTFilterOrFunc; kwargs...) =
    x -> iradon_threaded(x, f; kwargs...)


# function iradon_af(
#     sinog::AFMatrix{T};
#     rows = nothing,
#     cols = nothing,
#     filter::Union{Nothing,Function} = ram_lak,
#     α = 360,
#     α₀ = 0,
#     interpolation = nothing,
# ) where {T<:Real}
#     nd, nϕ = size(sinog)
#     if isnothing(rows)
#         rows = round(Int, nd / sqrt(2))
#     end
#     if isnothing(cols)
#         cols = rows
#     end
#     if isnothing(interpolation)
#         interpolation = mat -> LinearInterpolation(axes(mat), mat)
#     end
#     filtered = sinog
#     if !isnothing(filter)
#         af_filter = AFArray(filter(T, nd, nϕ))
#         filtered = fft(filtered, 1) .* af_filter
#         filtered = ifft(filtered, 1) .|> real
#     end
#     interp = interpolation(filtered)
#     image = AFArray(zeros(T, rows, cols))
#     m = length(image)
#     xs = AFArray(zeros(T, m))
#     ys = AFArray(zeros(T, m))
#     x₀::T = T(cols - 1) / 2
#     y₀::T = T(rows - 1) / 2
#     t₀::T = T(nd + 1) / 2
#
#     for k = 1:m
#         x::T = (k - 1) ÷ rows - x₀
#         y::T = (k - 1) % rows - y₀
#         xs[k] = x
#         ys[k] = y
#     end
#
#     Δϕ::T = deg2rad(α) / nϕ
#     ϕ₀::T = deg2rad(α₀)
#
#     @info "Computing inverse Radon transform..."
#
#     for iϕ = 1:nϕ
#         ϕ::T = ((iϕ - 1) * Δϕ + ϕ₀)
#         cϕ = cos(ϕ)
#         sϕ = sin(ϕ)
#         for k = 1:m
#             x = xs[k]
#             y = ys[k]
#             # To be consistent with our conventions should be '+'.
#             t = x * cϕ + y * sϕ + t₀
#             if 1 <= t <= nd
#                 image[k] += interp(t, iϕ)
#             end
#         end
#     end
#
#     image
# end
#
#
# function iradon_af_alt(
#     #sinog::AFMatrix{T};
#     sinog::Matrix{T};
#     rows = nothing,
#     cols = nothing,
#     filter::Union{Nothing,Function} = ram_lak,
#     α = 360,
#     α₀ = 0,
#     interpolation = nothing,
# ) where {T<:Real}
#     nd, nϕ = size(sinog)
#     if isnothing(rows)
#         rows = round(Int, nd / sqrt(2))
#     end
#     if isnothing(cols)
#         cols = rows
#     end
#     if isnothing(interpolation)
#         interpolation = mat -> LinearInterpolation(axes(mat), mat)
#     end
#     filtered = sinog
#     if !isnothing(filter)
#         #af_filter = AFArray(filter(T, nd, nϕ))
#         af_filter = filter(T, nd, nϕ)
#         filtered = fft(filtered, 1) .* af_filter
#         filtered = ifft(filtered, 1) .|> real
#     end
#     interp = interpolation(filtered)
#     #image = AFArray(zeros(T, rows, cols))
#     # cph = AFArray(Vector{T}(undef, nϕ))
#     # sph = AFArray(Vector{T}(undef, nϕ))
#     #m = length(image)
#     #m = rows * cols;
#     Δϕ::T = deg2rad(α) / nϕ
#     ϕ₀::T = deg2rad(α₀)
#     x₀::T = T(cols - 1) / 2
#     y₀::T = T(rows - 1) / 2
#     t₀::T = T(nd + 1) / 2
#
#     ϕs = T.(0:(nϕ-1)) * Δϕ .+ ϕ₀
#
#     tfs = zeros(T, 3, 2, nϕ)
#
#     for (iϕ, ϕ) in enumerate(ϕs)
#         c, s = cos(ϕ), sin(ϕ)
#         t₀ϕ = t₀ - x₀ * c - y₀ * s
#         tfs[:,:,iϕ] .= [[c zero(T)]; [s zero(T)]; [t₀ϕ iϕ]]
#     end
#
#     #tfsa = AFArray(tfs)
#     tfsa = tfs
#
#     @info "Computing inverse Radon transform..."
#
#     sum(
#         map(1:nϕ) do iϕ
#             tf = tfsa[:,:,iϕ]
#             transform(filtered, tf, rows, cols, AF_INTERP_BILINEAR, true)
#         end
#     )
# end


@inline function iradon(
    sinog::AbstractMatrix,
    g::AbstractParallelBeamGeometry,
    f::Optional{F} = nothing;
    kwargs...
) where {F <: CTFilterOrFunc}
    iradon_fast_threaded(sinog, g, f; kwargs...)
end


@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractFanBeamGeometry,
    f::Optional{F} = nothing;
    kwargs...
) where {F <: CTFilterOrFunc}
    g′, par_sinog = fan2para(sinog, geometry)
    par_sinog ↣ iradon(g′, f; kwargs...)
end


iradon(sinog::AbstractMatrix, g::AbstractGeometry; kwargs...) =
    iradon_fast_threaded(sinog, g; kwargs...)

iradon(sinog::AbstractMatrix, f::CTFilterOrFunc; kwargs...) =
    iradon_fast_threaded(sinog, f; kwargs...)

iradon(sinog::AbstractMatrix; kwargs...) = iradon_fast_threaded(sinog; kwargs...)

# iradon(sinog::AbstractMatrix; kwargs...) = iradon_threaded(sinog; kwargs...)
# iradon(sinog::AFMatrix; kwargs...) = iradon_af_alt(sinog; kwargs...)
iradon(; kwargs...) = x -> iradon(x; kwargs...)
iradon(g::AbstractGeometry, f::CTFilterOrFunc; kwargs...) = x -> iradon(x, g, f; kwargs...)
iradon(g::AbstractGeometry; kwargs...) = x -> iradon(x, g; kwargs...)
iradon(f::CTFilterOrFunc; kwargs...) = x -> iradon(x, f; kwargs...)


struct FBP{
    Geometry <: AbstractGeometry,
    Filter <: AbstractCTFilter
} <: AbstractIRadonAlgorithm
    geometry::Geometry
    filter::Filter
end


datatype(fbp::FBP) = datatype(fbp.geometry)
alg_geometry(fbp::FBP) = fbp.geometry
alg_params(::FBP) = nothing


FBP(g::AbstractGeometry) = FBP(g, RamLak())
FBP(g::AbstractGeometry, f::Function) = FBP(g, CTFilter(f))


iradon(data::AbstractMatrix, fbp::FBP; kwargs...) =
    iradon(data, fbp.geometry, fbp.filter; kwargs...)

(fbp::FBP)(data::AbstractMatrix; kwargs...) = iradon(data, fbp; kwargs...)
(fbp::FBP)(; kwargs...) = iradon(fbp; kwargs...)


reconstruct_image(sinog::CTSinogram, fbp::FBP; kwargs...) =
    sinog ↣ iradon(fbp; kwargs...)
