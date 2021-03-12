function fbpa(
    sinog::AbstractMatrix{T};
    rows::Optional{Integer} = nothing,
    cols::Optional{Integer} = nothing,
    filter::Optional{F} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    ν::Real = 1,
    background::Optional{<:Real} = nothing,
    interpolation::Optional{Interp} = nothing,
    progress::Bool = true,
) where {
    T <: Real,
    F <: Filters.CTFilterOrFunc,
    Interp <: Union{Function,AbstractInterp2DOrNone}
}
    M = typeof(sinog)
    nd, nϕ = size(sinog)
    rows = maybe(cols, rows)
    rows = maybe(round(Int, nd / √2), rows)
    cols = maybe(rows, cols)
    filtered = apply(maybe(Filters.RamLak(), filter)) do f
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
    p = Progress(
        length(image);
        dt=0.2,
        desc="Computing inverse Radon transform...",
        enabled=progress,
    )
    Threads.@threads for k ∈ eachindex(xys)
        x, y = xys[k]
        @inbounds image[k] = sum(eachindex(scϕs)) do iϕ
            sϕ, cϕ = scϕs[iϕ]
            # To be consistent with our conventions should be '+'.
            t = x * cϕ + y * sϕ + t₀
            t ∈ 1..nd ? interp(t, T(iϕ)) : z
        end
        next!(p)
    end
    image |> CTTomogram{M}
end
