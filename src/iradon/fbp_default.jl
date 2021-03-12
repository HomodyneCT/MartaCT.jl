function fbp(
    sinog::AbstractMatrix{T};
    rows::Optional{Integer} = nothing,
    cols::Optional{Integer} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    ν::Real = 1,
    filter::CTFilter = RamLak(),
    background::Real = zero(T),
    interpolation::Interp = BilinearInterpolation(),
    progress::Bool = true,
) where {T <: Real,Interp <: AbstractInterp2DOrNone}
    nd, nϕ = size(sinog)
    rows = maybe(cols, rows)
    rows = maybe(round(Int, nd / √2), rows)
    cols = maybe(rows, cols)
    filtered = apply(filter) do f
        filter_freq = fft(sinog, 1) .* f(T, nd, nϕ)
        ifft(filter_freq, 1) |> real
    end
    interp = interpolation(filtered)
    Δϕ::T = deg2rad(α) / nϕ
    ϕ₀::T = deg2rad(α₀)
    t₀::T = T(nd + 1) / 2
    sθ, cθ = sincos(atan(T(rows), T(cols)))
    x₀, y₀ = t₀ * cθ / ν, t₀ * sθ / ν
    scϕs = map(sincos, range(ϕ₀; step = Δϕ, length = nϕ))
    xs = linspace(-x₀..x₀, cols)
    ys = linspace(-y₀..y₀, rows)
    xys = Vector{NTuple{2,T}}(undef, rows * cols)
    @inbounds @simd for k ∈ eachindex(xys)
        xys[k] = xs[(k - 1) ÷ rows + 1], ys[(k - 1) % rows + 1]
    end
    image = similar(sinog, rows, cols)
    fill!(image, T(background))
    temp_images = fill(image, 1)
    Threads.resize_nthreads!(temp_images)
    p = @iradonprogress nϕ progress
    Threads.@threads for iϕ ∈ eachindex(scϕs)
        sϕ, cϕ = scϕs[iϕ]
        id = Threads.threadid()
        @inbounds img = temp_images[id]
        @inbounds @simd for k ∈ eachindex(xys)
            x, y = xys[k]
            # To be consistent with our conventions should be '+'.
            t::T = x * cϕ + y * sϕ + t₀
            if t ∈ 1..nd
                img[k] += interp(t, T(iϕ))
            end
        end
        next!(p)
    end
    CTTomogram(sum!(image, temp_images))
end


function fbp(
    sinog::AbstractMatrix,
    geometry::AbstractParallelBeamGeometry,
    filter::CTFilter = RamLak();
    kwargs...
)
    rows = geometry.rows
    cols = geometry.cols
    α = geometry.α
    α₀ = geometry.α₀
    fbp(sinog; rows, cols, α, α₀, filter, kwargs...)
end








function iradon_threaded(
    sinog::AbstractMatrix{T},
    geometry::AbstractParallelBeamGeometry,
    filter::Optional{F} = nothing;
    kwargs...
) where {T <: Real,F <: Filters.CTFilterOrFunc}
    rows = geometry.rows
    cols = geometry.cols
    α = geometry.α
    α₀ = geometry.α₀
    iradon_threaded(sinog; rows, cols, α, α₀, filter, kwargs...)
end


function iradon_threaded(
    sinog::AbstractMatrix,
    filter::Filters.CTFilterOrFunc;
    geometry::Optional{G} = nothing,
    kwargs...
) where {G <: AbstractParallelBeamGeometry}
    isnothing(geometry) && return iradon_threaded(sinog; filter, kwargs...)
    iradon_threaded(sinog, geometry, filter; kwargs...)
end


iradon(sinog::AbstractMatrix; kwargs...) = default_iradon()(sinog; kwargs...)


@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractParallelBeamGeometry,
    filter::Optional{F} = nothing;
    kwargs...
) where {F <: Filters.CTFilterOrFunc}
    rows = geometry.rows
    cols = geometry.cols
    α = geometry.α
    α₀ = geometry.α₀
    iradon(sinog; rows, cols, α, α₀, filter, kwargs...)
end


@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractFanBeamGeometry,
    f::Optional{F} = nothing;
    kwargs...
) where {F <: Filters.CTFilterOrFunc}
    g′, par_sinog = fan2para(sinog, geometry)
    par_sinog ↣ iradon(g′, f; kwargs...)
end


iradon(sinog::AbstractMatrix, g::AbstractGeometry; kwargs...) =
    default_iradon()(sinog, g; kwargs...)

iradon(sinog::AbstractMatrix, f::Filters.CTFilterOrFunc; kwargs...) =
    default_iradon()(sinog, f; kwargs...)

iradon(; kwargs...) = x -> iradon(x; kwargs...)
iradon(g::AbstractGeometry, f::Filters.CTFilterOrFunc; kwargs...) =
    x -> iradon(x, g, f; kwargs...)
iradon(g::AbstractGeometry; kwargs...) = x -> iradon(x, g; kwargs...)
iradon(f::Filters.CTFilterOrFunc; kwargs...) = x -> iradon(x, f; kwargs...)


struct FBP{F<:AbstractCTFilter} <: AbstractIRadonAlgorithm
    filter::F
end

FBP() = FBP(RamLak())
FBP(f::Function) = FBP(CTFilter(f))


iradon(data::AbstractMatrix, fbp::FBP; kwargs...) =
    iradon(data, fbp.geometry, fbp.filter; kwargs...)

(fbp::FBP)(data::AbstractMatrix; kwargs...) = iradon(data, fbp; kwargs...)
(fbp::FBP)(; kwargs...) = iradon(fbp; kwargs...)
