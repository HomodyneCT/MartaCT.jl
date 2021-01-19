module FanBeam

export para2fan, fan2para

using ..Monads, ..Applicative
using ..CTImages: CTSinogram
using ..Geometry
using ..Interpolation: interpolate, AbstractBilinearInterpolation


"""
    para2fan(
        sinog_para::AbstractMatrix{T},
        fbg::FanBeamGeometry{T};
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
- `interpolation`: interpolation strategy; it should be a function
    taking a matrix as input and returning a function of the indices to
    get the interpolated value. By default it is a bilinear interpolation.

See Also: [`fan2para`](@ref)
"""
function para2fan(
    sinog_para::AbstractMatrix{T},
    fbg::FanBeamGeometry{T};
    interpolation::Optional{Interp} = nothing,
) where {T <: Real, Interp <: Union{Function,AbstractBilinearInterpolation}}
    nd, nϕ = num_det(fbg), num_proj(fbg)

    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(sinog_para)
    sinog_fan = zeros(T, nd, nϕ)

    D::T = f2iso(fbg)
    Δβ::T = (deg2rad ∘ T ∘ scan_angle)(fbg) / T(nϕ-1) # This should be 2π / (nβ-1)?
    β₀::T = (deg2rad ∘ T ∘ start_angle)(fbg) + T(1)
    Δϕ::T = T(2π) / T(nϕ-1) # This should be α / (nϕ-1)?
    γ::T = (T ∘ fan_angle)(fbg)
    Δγ::T = γ / T(nd)
    x′max::T = T(D) * sin(γ/2)
    Δx′::T = 2x′max / T(nd)
    center = center_channel(fbg)
    γ₀::T = center + 1
    x′₀::T = center + 1

    function compute_value(x′::T, ϕ::T)
        if ϕ < 0
            ϕ += π
            x′ = -x′
        end
        if ϕ > 2π
            ϕ -= π
            x′ = -x′
        end
        x′ = x′ / Δx′ + x′₀
        ϕ = mod2pi(ϕ) / Δϕ + 1 # need +1 in order to be in 1:nϕ
        if 1 <= x′ <= nd && 1 <= ϕ <= nϕ
            return interp(x′, ϕ)
        end
        return 0
    end

    Threads.@threads for iγ = 1:nd
        γ′::T = (T(iγ) - γ₀) * Δγ
        x′::T = T(D) * sin(γ′)
        for iβ = 1:nϕ
            β::T = T(iβ - β₀) * Δβ
            ϕ::T = β - γ′
            @inbounds sinog_fan[iγ, iβ] = compute_value(x′, ϕ)
        end
    end

    fbg, CTSinogram(sinog_fan)
end



function para2fan(
    sinog_para::AbstractMatrix{T},
    pbg::ParallelBeamGeometry{T};
    D,
    D′ = nothing,
    γ = nothing,
    δ = 1,
    kwargs...,
) where {T<:Real}
    @assert (pbg.nd, pbg.nϕ) == size(sinog_para) "Sinogram size $(size(sinog_para)) should match geometry ($(pbg.nd),$(pbg.nϕ))"

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

function para2fan(img::CTSinogram, g::AbstractGeometry; kwargs...)
    mbind(img) do x
        g′, res = para2fan(x, g; kwargs...)
        g′, res
    end
end


"""
    fan2para(
        sinog_fan::AbstractMatrix{T},
        fbg::FanBeamGeometry{T};
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
- `interpolation`: interpolation strategy; it should be a function
    taking a matrix as input and returning a function of the indices to
    get the interpolated value. By default it is a bilinear interpolation.

See Also: [`para2fan`](@ref)
"""
function fan2para(
    sinog_fan::AbstractMatrix{T},
    fbg::FanBeamGeometry{T,DefaultTomograph};
    interpolation::Optional{Interp} = nothing,
) where {T <: Real, Interp <: Union{Function,AbstractBilinearInterpolation}}
    nd, nϕ = fbg.nd, fbg.nϕ

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

    interpolation = maybe(interpolate, interpolation)
    interp = interpolation(sinog_fan)
    sinog_para = zeros(T, nd, nϕ)
    D::T = f2iso(fbg)
    γ::T = fan_angle(fbg)
    Δϕ::T = (deg2rad ∘ T ∘ scan_angle)(fbg) / T(nϕ - 1)
    Δβ::T = 2π / T(nϕ - 1)
    β₀::T = (deg2rad ∘ T ∘ start_angle)(fbg)
    Δγ::T = γ / T(nd)
    x′max::T = D * sin(γ/2)
    Δx′::T = 2x′max / T(nd)
    center = center_channel(pbg)
    x′₀::T = center + 1 # support for custom center channel
    @assert (x′₀ * Δx′ / D) <= 1 "Fan beam parameters incompatible |x′₀/D| = $(x′₀ * Δx′ / D) which should be less than 1"
    γ₀::T = center + 1 # support for custom center channel.

    function compute_value(γ′::T, β::T)
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
        γ′ = γ′ / Δγ + γ₀
        β = mod2pi(β) / Δβ + 1 # need +1 in order to be in 1:nβ
        if 1 <= γ′ <= nd && 1 <= β <= nϕ
            return interp(γ′, β)
        end
        return 0
    end

    Threads.@threads for ix = 1:nd
        x′::T = (ix - x′₀) * Δx′
        γ′::T = asin(x′ / D)
        for iϕ = 1:nϕ
            ϕ::T = (iϕ - 1) * Δϕ
            β::T = ϕ - γ′ + β₀ # Default tomograph geometry requires '-'
            @inbounds sinog_para[ix, iϕ] = compute_value(γ′, β)
        end
    end

    pbg, CTSinogram(sinog_para)
end


fan2para(g::AbstractGeometry; kwargs...) = x -> fan2para(x, g; kwargs...)

function fan2para(img::CTSinogram, g::AbstractGeometry; kwargs...)
    mbind(img) do x
        g′, res = fan2para(x, g; kwargs...)
        g′, res
    end
end

end # module
