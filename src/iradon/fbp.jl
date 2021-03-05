function iradon_fast_threaded(
    sinog::AbstractMatrix{T};
    rows::Optional{Integer} = nothing,
    cols::Optional{Integer} = nothing,
    filter::Optional{F} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    ν::Real = 1,
    background::Optional{<:Real} = nothing,
    interpolation::Optional{Interp} = nothing,
) where {
    T <: Real,
    F <: CTFilterOrFunc,
    Interp <: Union{Function,AbstractInterp2DOrNone}
}
    M = typeof(sinog)
    nd, nϕ = size(sinog)
    rows = maybe(cols, rows)
    rows = maybe(round(Int, nd / √2), rows)
    cols = maybe(rows, cols)
    filtered = apply(maybe(RamLak(), filter)) do f
        filter_freq = fft(sinog, 1) .* f(T, nd, nϕ)
        ifft(filter_freq, 1) |> real
    end
    interpolation = maybe(interpolate, interpolation)
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
    z::T = maybe(zero(T), background)
    temp_images = fill(fill(z, rows, cols), 1)
    Threads.resize_nthreads!(temp_images)
    @info "Computing inverse Radon transform..."
    p = Progress(nϕ, 0.2)
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
    temp_images |> sum |> CTTomogram{M}
end


function iradon_fast_threaded_alt(
    sinog::AbstractMatrix{T};
    rows::Optional{Integer} = nothing,
    cols::Optional{Integer} = nothing,
    filter::Optional{F} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    ν::Real = 1,
    background::Optional{<:Real} = nothing,
    interpolation::Optional{Interp} = nothing,
) where {
    T <: Real,
    F <: CTFilterOrFunc,
    Interp <: Union{Function,AbstractInterp2DOrNone}
}
    M = typeof(sinog)
    nd, nϕ = size(sinog)
    rows = maybe(cols, rows)
    rows = maybe(round(Int, nd), rows)
    cols = maybe(rows, cols)
    l = min(rows, cols)
    filtered = apply(maybe(RamLak(), filter)) do f
        filter_freq = fft(sinog, 1) .* f(T, nd, nϕ)
        ifft(filter_freq, 1) |> real
    end
    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(filtered)
    Δϕ::T = deg2rad(α) / nϕ
    ϕ₀::T = deg2rad(α₀)
    t₀::T = T(nd + 1) / 2
    x₀, y₀ = (cols - l) ÷ 2, (rows - l) ÷ 2
    scϕs = map(sincos, range(ϕ₀; step = Δϕ, length = nϕ))
    xs = linspace(-t₀..t₀, l) / ν
    ys = linspace(-t₀..t₀, l) / ν
    indices = Vector{NTuple{2,Int}}(undef, l^2)
    xys = Vector{NTuple{2,T}}(undef, l^2)
    @inbounds @simd for k ∈ eachindex(xys)
        ix, iy = (k - 1) ÷ l + 1, (k - 1) % l + 1
        indices[k] = x₀ + ix, y₀ + iy
        xys[k] = xs[ix], ys[iy]
    end
    z::T = maybe(zero(T), background)
    temp_images = fill(fill(z, rows, cols), 1)
    Threads.resize_nthreads!(temp_images)
    @info "Computing inverse Radon transform..."
    p = Progress(nϕ, 0.2)
    Threads.@threads for iϕ ∈ eachindex(scϕs)
        sϕ, cϕ = scϕs[iϕ]
        id = Threads.threadid()
        @inbounds img = temp_images[id]
        @inbounds @simd for k ∈ eachindex(xys)
            ix, iy = indices[k]
            x, y = xys[k]
            # To be consistent with our conventions should be '+'.
            t::T = x * cϕ + y * sϕ + t₀
            if t ∈ 1..nd
                img[iy,ix] += interp(t, T(iϕ))
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
) where {T <: Real,F <: CTFilterOrFunc}
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
    isnothing(geometry) && return iradon_fast_threaded(sinog; filter, kwargs...)
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
    ν::Real = 1,
    background::Optional{<:Real} = nothing,
    interpolation::Optional{Interp} = nothing,
) where {
    T <: Real,
    F <: CTFilterOrFunc,
    Interp <: Union{Function,AbstractInterp2DOrNone}
}
    M = typeof(sinog)
    nd, nϕ = size(sinog)
    rows = maybe(cols, rows)
    rows = maybe(round(Int, nd / √2), rows)
    cols = maybe(rows, cols)
    filtered = apply(maybe(RamLak(), filter)) do f
        filter_freq = fft(sinog, 1) .* f(T, nd, nϕ)
        ifft(filter_freq, 1) |> real
    end
    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(filtered)
    Δϕ::T = deg2rad(α) / nϕ
    ϕ₀::T = deg2rad(α₀)
    t₀::T = T(nd + 1) / 2
    sθ, cθ = sincos(atan(T(rows), T(cols)))
    x₀, y₀ = t₀ * cθ / ν, t₀ * sθ / ν
    xs = linspace(-x₀..x₀, cols)
    ys = linspace(-y₀..y₀, rows)
    xys = Vector{NTuple{2,T}}(undef, rows * cols)
    @inbounds @simd for k ∈ eachindex(xys)
        xys[k] = xs[(k - 1) ÷ rows + 1], ys[(k - 1) % rows + 1]
    end
    scϕs = map(sincos, range(ϕ₀; step = Δϕ, length = nϕ))
    z::T = maybe(zero(T, background))
    image = fill(z, rows, cols)
    @info "Computing inverse Radon transform..."
    p = Progress(length(image), 0.2)
    Threads.@threads for k ∈ eachindex(xys)
        x, y = xys[k]
        @inbounds image[k] = sum(eachindex(scϕs)) do iϕ
            sϕ, cϕ = scϕs[iϕ]
            # To be consistent with our conventions should be '+'.
            t::T = x * cϕ + y * sϕ + t₀
            t ∈ 1..nd ? interp(t, T(iϕ)) : z
        end
        next!(p)
    end
    image |> CTTomogram{M}
end


function iradon_threaded_alt(
    sinog::AbstractMatrix{T};
    rows::Optional{Integer} = nothing,
    cols::Optional{Integer} = nothing,
    filter::Optional{F} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    ν::Real = 1,
    background::Optional{<:Real} = nothing,
    interpolation::Optional{Interp} = nothing,
) where {
    T <: Real,
    F <: CTFilterOrFunc,
    Interp <: Union{Function,AbstractInterp2DOrNone}
}
    M = typeof(sinog)
    nd, nϕ = size(sinog)
    rows = maybe(cols, rows)
    rows = maybe(round(Int, nd), rows)
    cols = maybe(rows, cols)
    l = min(rows, cols)
    filtered = apply(maybe(RamLak(), filter)) do f
        filter_freq = fft(sinog, 1) .* f(T, nd, nϕ)
        ifft(filter_freq, 1) |> real
    end
    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(filtered)
    Δϕ::T = deg2rad(α) / nϕ
    ϕ₀::T = deg2rad(α₀)
    t₀::T = T(nd + 1) / 2
    x₀, y₀ = (cols - l) ÷ 2, (rows - l) ÷ 2
    xs = linspace(-t₀..t₀, l) / ν
    ys = linspace(-t₀..t₀, l) / ν
    xys = Vector{NTuple{2,T}}(undef, l^2)
    indices = Vector{NTuple{2,Int}}(undef, l^2)
    @inbounds @simd for k ∈ eachindex(xys)
        ix, iy = (k - 1) ÷ l + 1, (k - 1) % l + 1
        indices[k] = x₀ + ix, y₀ + iy
        xys[k] = xs[ix], ys[iy]
    end
    scϕs = map(sincos, range(ϕ₀; step = Δϕ, length = nϕ))
    z::T = maybe(zero(T), background)
    image = fill(z, rows, cols)
    @info "Computing inverse Radon transform..."
    p = Progress(length(xys), 0.2)
    Threads.@threads for k ∈ eachindex(xys)
        ix, iy = indices[k]
        x, y = xys[k]
        @inbounds image[iy, ix] = sum(eachindex(scϕs)) do iϕ
            sϕ, cϕ = scϕs[iϕ]
            # To be consistent with our conventions should be '+'.
            t::T = x * cϕ + y * sϕ + t₀
            t ∈ 1..nd ? interp(t, T(iϕ)) : z
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
) where {T <: Real,F <: CTFilterOrFunc}
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


const _default_iradon_ref = Ref{Function}(iradon_fast_threaded)
default_iradon() = _default_iradon_ref[]
default_iradon(other::Function) = _default_iradon_ref[] = other
default_iradon(other::Symbol) = @eval _default_iradon_ref[] = $other


iradon(sinog::AbstractMatrix; kwargs...) = default_iradon()(sinog; kwargs...)


@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractParallelBeamGeometry,
    filter::Optional{F} = nothing;
    kwargs...
) where {F <: CTFilterOrFunc}
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
) where {F <: CTFilterOrFunc}
    g′, par_sinog = fan2para(sinog, geometry)
    par_sinog ↣ iradon(g′, f; kwargs...)
end


iradon(sinog::AbstractMatrix, g::AbstractGeometry; kwargs...) =
    default_iradon()(sinog, g; kwargs...)

iradon(sinog::AbstractMatrix, f::CTFilterOrFunc; kwargs...) =
    default_iradon()(sinog, f; kwargs...)

iradon(; kwargs...) = x -> iradon(x; kwargs...)
iradon(g::AbstractGeometry, f::CTFilterOrFunc; kwargs...) = x -> iradon(x, g, f; kwargs...)
iradon(g::AbstractGeometry; kwargs...) = x -> iradon(x, g; kwargs...)
iradon(f::CTFilterOrFunc; kwargs...) = x -> iradon(x, f; kwargs...)


struct FBP{Geometry <: AbstractGeometry,Filter <: AbstractCTFilter} <: AbstractIRadonAlgorithm
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
