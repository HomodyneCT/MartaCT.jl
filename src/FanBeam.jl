module FanBeam

export para2fan, fan2para

using ..Monads
using ..CTImages: CTSinogram
using ..Geometry
using ..Interpolation: interpolate, AbstractInterp2DOrNone
using ..Utils: linspace, ORI, _wrap_angle
using IntervalSets


"""
    para2fan(
        sinog_para::AbstractMatrix{T},
        fbg::FanBeamGeometry;
        <keyword arguments>
    ) where {T<:Real} -> fbg, sinog_fan

Convert given sinogram `sinog_para` from parallel beam geometry
to fan beam projections.

The function returns a tuple where `fbg` is the new geometry
with fan beam projections parameters and `sinog_fan` is the
converted sinogram.

# Arguments
- `sinog_para`: the input sinogram.
- `pbg`: the geometry parameters.
- `D`: focal spot to ISO distance.
- `D′=nothing`: focal spot to detectors distance; if not specified
    and also `γ` is not specified, is defined so that the total
    fan angle is 1 radian.
- `γ=nothing`: total fan angle in degrees.
- `δ=1`: detectors spacing (cell size).
- `background=nothing`: background value to be used, defaults to 0.
- `interpolation`: interpolation strategy; it should be a function
    taking a matrix as input and returning a function of the indices to
    get the interpolated value. By default it is a bilinear interpolation.

See Also: [`fan2para`](@ref)
"""
function para2fan(
    sinog_para::AbstractMatrix{T},
    fbg::FanBeamGeometry;
    background::Optional{Real} = nothing,
    interpolation::Optional{Interp} = nothing,
) where {T <: Real, Interp <: AbstractInterp2DOrNone}
    nd, nϕ = num_det(fbg), num_proj(fbg)

    @assert (nd, nϕ) == size(sinog_para) "Sinogram size $(size(sinog_para)) should match geometry ($nd,$nϕ)"

    D::T = f2iso(fbg)
    Δβ::T = deg2rad(scan_angle(fbg)) # This should be 2π?
    βs = linspace(T, ORI(0..Δβ), nϕ)
    β₀::T = deg2rad(start_angle(fbg)) + 1
    δϕi::T = nϕ / 2π # This should be α / nϕ?
    γ::T = fan_angle(fbg)
    δγ::T = γ / nd
    x′max::T = D * sin(γ/2)
    δx′i::T = nd / 2x′max
    center = center_channel(fbg)
    γ₀::T = center * δγ
    γs = linspace(T, ORI(-γ₀..γ₀), nd)
    x′₀::T = center + 1
    nthalf = nϕ + 1//2

    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(sinog_para)
    z::T = maybe(zero(T), background)
    sinog_fan = similar(sinog_para, nd, nϕ)
    fill!(sinog_fan, z)

    @inline function compute_value(x′::T, ϕ::T)
        if ϕ < 0
            ϕ += π
            x′ = -x′
        end
        if ϕ > 2π
            ϕ -= π
            x′ = -x′
        end
        x′ = x′ * δx′i + x′₀
        ϕ = mod2pi(ϕ) * δϕi + 1 # need +1 in order to be in 1:nϕ
        return x′ ∈ 1..nd && ϕ ∈ 1..nϕ+1 ?
            interp(x′, _wrap_angle(ϕ, nϕ, nthalf)) : z
    end

    Threads.@threads for iγ ∈ 1:nd
        @inbounds γ′ = γs[iγ]
        x′ = D * sin(γ′)
        @inbounds @simd for iβ ∈ 1:nϕ
            β = βs[iβ]
            ϕ = β - γ′
            sinog_fan[iγ, iβ] = compute_value(x′, ϕ)
        end
    end

    fbg, CTSinogram(sinog_fan)
end


function para2fan(
    sinog_para::AbstractMatrix{T},
    pbg::ParallelBeamGeometry;
    D::Real = 500,
    D′::Optional{Real} = nothing,
    D1::Optional{Real} = nothing,
    γ::Optional{Real} = nothing,
    gamma::Optional{Real} = nothing,
    δ::Optional{Real} = one(T),
    dx::Optional{Real} = one(T),
    kwargs...,
) where {T<:Real}
    D′ = maybe(D′, D1)
    γ = maybe(γ, gamma)
    δ = maybe(δ, dx)
    fbg = FanBeamGeometry(
        T,
        pbg.ct;
        pbg.nϕ,
        pbg.nd,
        pbg.rows,
        pbg.cols,
        D,
        D′,
        γ,
        δ,
        pbg.α,
        pbg.α₀,
        pbg.center,
    )
    para2fan(sinog_para, fbg; kwargs...)
end


para2fan(g::AbstractGeometry; kwargs...) = x -> para2fan(x, g; kwargs...)


"""
    fan2para(
        sinog_fan::AbstractMatrix{T},
        fbg::FanBeamGeometry;
        <keyword arguments>
    ) where {T<:Real} -> pbg, sinog_para

Convert given sinogram `sinog_fan` from fan beam geometry
to parallel beam projections.

It is assumed that the detectors are arranged on a circle arc,
not flat.

The function returns a tuple where `fbg` is the new geometry
with fan beam projections parameters and `sinog_fan` is the
converted sinogram.

# Arguments
- `sinog_para`: the input sinogram.
- `pbg`: the parallel geometry parameters.
- `background=nothing`: background value to be used, defaults to 0.
- `interpolation`: interpolation strategy; it should be a function
    taking a matrix as input and returning a function of the indices to
    get the interpolated value. By default it is a bilinear interpolation.

See Also: [`para2fan`](@ref)
"""
function fan2para(
    sinog_fan::AbstractMatrix{T},
    fbg::FanBeamGeometry{T,DefaultTomograph};
    background::Optional{Real} = nothing,
    interpolation::Optional{Interp} = nothing,
) where {T <: Real, Interp <: AbstractInterp2DOrNone}
    nd, nϕ = num_det(fbg), num_proj(fbg)

    @assert (nd, nϕ) == size(sinog_fan) "Sinogram size $(size(sinog_fan)) should match geometry ($nd,$nϕ)"

    pbg = ParallelBeamGeometry(
        T,
        fbg.ct;
        nϕ,
        nd,
        fbg.rows,
        fbg.cols,
        fbg.α,
        fbg.center,
    )

    D::T = f2iso(fbg)
    γ::T = fan_angle(fbg)
    Δϕ::T = deg2rad(scan_angle(fbg))
    ϕs = linspace(T, ORI(0..Δϕ), nϕ)
    δβi::T = nϕ / 2π
    β₀::T = deg2rad(start_angle(fbg))
    δγi::T = nd / γ
    x′max::T = D * sin(γ/2)
    δx′::T = 2x′max / nd
    center::T = center_channel(pbg)
    x′₀::T = center * δx′
    x′₁::T = (nd - center) * δx′
    xs = linspace(T, -x′₀..x′₁, nd) # support for custom center channel
    @assert (x′₀ / D) <= 1 "Fan beam parameters incompatible |x′₀/D| = $(x′₀ / D) which should be less than 1"
    γ₀::T = center + 1 # support for custom center channel.

    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(sinog_fan)
    z::T = maybe(zero(T), background)
    sinog_para = similar(sinog_fan, nd, nϕ)
    fill!(sinog_para, z)
    nthalf = nϕ + 1//2

    @inline function compute_value(γ′::T, β::T)
        if β < 0
            # β += 2γ′ + π
            # γ′ = -γ′
            β += π - 2γ′ # We use by default the symmetry at π instead of 2π
            γ′ = -γ′
        end
        if β > 2π
            # β += 2γ′ - π
            # γ′ = -γ′
            β -= π + 2γ′
            γ′ = -γ′
        end
        γ′′ = γ′ * δγi + γ₀
        β = mod2pi(β) * δβi + 1 # need +1 in order to be in 1:nβ
        return γ′′ ∈ 1..nd && β ∈ 1..nϕ+1 ?
            interp(γ′′, _wrap_angle(β, nϕ, nthalf)) : z
    end

    Threads.@threads for ix ∈ 1:nd
        @inbounds x′ = xs[ix]
        γ′ = asin(x′ / D)
        @inbounds @simd for iϕ ∈ 1:nϕ
            ϕ = ϕs[iϕ]
            #β::T = ϕ - γ′ + β₀ # Default tomograph geometry requires '-'
            β = ϕ + γ′ + β₀ # It does not!
            sinog_para[ix, iϕ] = compute_value(γ′, β)
        end
    end

    pbg, CTSinogram(sinog_para)
end


fan2para(g::AbstractGeometry; kwargs...) = x -> fan2para(x, g; kwargs...)

end # module
