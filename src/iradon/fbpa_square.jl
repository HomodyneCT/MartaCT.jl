struct FBPAFFTSquare{F<:AbstractCTFilter} <: AbstractFBP
    filter::F
end

FBPAFFTSquare() = FBPAFFTSquare(RamLak())
FBPAFFTSquare(f::Function) = FBPAFFTSquare(CTFilter(f))

@defiradonalgfn FBPAFFTSquare fbpa_fft_square


function fbpa_fft_square(
    sinog::AbstractMatrix{T},
    xs::AbstractVector{U1},
    ys::AbstractVector{U2},
    ϕs::ClosedInterval = 0..2π;
    background::Optional{U3} = nothing,
    filter::Optional{F} = nothing,
    interpolation::Optional{Interp} = nothing,
    progress::Bool = true,
) where {
    T <: Real,
    U1 <: Real,
    U2 <: Real,
    U3 <: Real,
    F <: AbstractCTFilter,
    Interp <: AbstractInterp2DOrNone,
}
    nd, nϕ = size(sinog)
    cols = length(xs)
    rows = length(ys)
    l = min(rows, cols)
    filtered = apply(maybe(Filters.RamLak(), filter)) do f
        filter_freq = fft(sinog, 1) .* f(T, nd, nϕ)
        ifft(filter_freq, 1) |> real
    end
    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(filtered)
    t₀::T = (nd + 1) / 2
    x₀, y₀ = (cols - l) ÷ 2, (rows - l) ÷ 2
    xys = Vector{NTuple{2,T}}(undef, l^2)
    indices = Vector{NTuple{2,Int}}(undef, l^2)
    @inbounds @simd for k ∈ eachindex(xys)
        ix, iy = (k - 1) ÷ l + 1, (k - 1) % l + 1
        indices[k] = x₀ + ix, y₀ + iy
        xys[k] = xs[ix], ys[iy]
    end
    scϕs = sincos.(linspace(T, ORI(ϕs), nϕ))
    z::T = maybe(zero(T), background)
    image = similar(sinog, rows, cols)
    fill!(image, z)
    p = @iradonprogress length(xys) progress
    Threads.@threads for k ∈ eachindex(xys)
        ix, iy = indices[k]
        x, y = xys[k]
        @inbounds image[iy, ix] = sum(eachindex(scϕs)) do iϕ
            sϕ, cϕ = scϕs[iϕ]
            # To be consistent with our conventions should be '+'.
            t = x * cϕ + y * sϕ + t₀
            t ∈ 1..nd ? interp(t, T(iϕ)) : z
        end
        next!(p)
    end
    CTTomogram(image)
end


@defsquarekwfn fbpa_fft_square
@defiradonfngeom fbpa_fft_square
