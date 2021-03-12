function fbpa_square(
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
    rows = maybe(round(Int, nd), rows)
    cols = maybe(rows, cols)
    l = min(rows, cols)
    filtered = apply(maybe(Filters.RamLak(), filter)) do f
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
    p = Progress(
        length(xys);
        dt=0.2,
        desc="Computing inverse Radon transform...",
        enabled=progress,
    )
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
