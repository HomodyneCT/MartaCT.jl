struct FBP{F<:AbstractCTFilter} <: AbstractFBP
    filter::F
end

FBP() = FBP(RamLak())
FBP(f::Function) = FBP(CTFilter(f))

@defiradonalgfn FBP fbp_fft


function fbp_fft(
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
    filtered = apply(maybe(RamLak(), filter)) do f
        filter_freq = fft(sinog, 1) .* f(T, nd, nϕ)
        ifft(filter_freq, 1) |> real
    end
    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(filtered)
    t₀::T = (nd + 1) / 2
    scϕs = sincos.(linspace(T, ORI(ϕs), nϕ))
    xys = Vector{NTuple{2,T}}(undef, rows * cols)
    @inbounds @simd for k ∈ eachindex(xys)
        xys[k] = xs[(k - 1) ÷ rows + 1], ys[(k - 1) % rows + 1]
    end
    z::T = maybe(zero(T), background)
    temp_images = fill(fill(z, rows, cols), 1)
    Threads.resize_nthreads!(temp_images)
    p = @iradonprogress nϕ progress
    Threads.@threads for iϕ ∈ eachindex(scϕs)
        sϕ, cϕ = scϕs[iϕ]
        id = Threads.threadid()
        @inbounds img = temp_images[id]
        @inbounds @simd for k ∈ eachindex(xys)
            x, y = xys[k]
            # To be consistent with our conventions should be '+'.
            t = x * cϕ + y * sϕ + t₀
            if t ∈ 1..nd
                img[k] += interp(t, iϕ)
            end
        end
        next!(p)
    end
    image = similar(sinog, rows, cols)
    CTTomogram(sum!(image, temp_images))
end

@defdiagkwfn fbp_fft
@defiradonfngeom fbp_fft
