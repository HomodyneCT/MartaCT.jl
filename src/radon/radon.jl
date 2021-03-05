"""radon_std(image::AbstractMatrix; <keyword arguments>)

Compute the Radon transform of `image` inside a circle of
radius `hypot(rows,cols)/2` where `rows` and `cols` are the
dimensions of `image`.

See Also: [`radon_alt`](@ref)
"""
function radon_std(
    image::AbstractMatrix{T};
    nd::Optional{Int} = nothing,
    nϕ::Optional{Int} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    ν::Real = 1,
    background::Optional{<:Real} = nothing,
    rescaled::Bool = true,
    interpolation::Optional{Interp} = nothing,
) where {T <: Real,Interp <: Union{Function,AbstractInterp2DOrNone}}
    M = typeof(image)
    rows, cols = size(image)
    nd = isnothing(nd) ? round(Int, hypot(rows, cols)) : nd
    nϕ = isnothing(nϕ) ?
        2 * (round(Int, (rows + 1) * (cols + 1) / nd) ÷ 2) + 1 : nϕ
    x′₀::T = hypot(rows, cols) / 2
    sθ, cθ = atan(rows, cols) |> sincos
    x₀::T = x′₀ * cθ + 1
    y₀::T = x′₀ * sθ + 1
    ϕ₀::T = deg2rad(α₀)
    Δϕ::T = deg2rad(α) / nϕ
    x′s = linspace(-x′₀..x′₀, nd) * ν
    scϕs = map(sincos, range(ϕ₀; step = Δϕ, length = nϕ))
    z::T = maybe(zero(T), background)
    rmat = fill(z, nd, nϕ)
    indices = similar(vec(rmat), NTuple{2,Int})
    x′ϕs = similar(indices, NTuple{2,T})
    foreach(eachindex(rmat)) do k
        ix = (k - 1) % nd
        iϕ = (k - 1) ÷ nd
        @inbounds x′::T = x′s[ix + 1]
        @inbounds s, c = scϕs[iϕ + 1]
        @inbounds indices[k] = k, iϕ
        @inbounds x′ϕs[k] = x′ * c, x′ * s
    end
    rimage = rescaled ? rescale(image) : image
    interp = isnothing(interpolation) ?
        interpolate(rimage) : interpolation(rimage)
    @info "Computing Radon transform..."
    p = Progress(length(rmat), 0.2)
    Threads.@threads for (k, iϕ) ∈ indices
        @inbounds x′x, x′y = x′ϕs[k]
        prex, prey = x′x + x₀, x′y + y₀
        o = iϕ * nd
        @inbounds rmat[k] = sum(view(x′ϕs, o+1:o+nd)) do (y′y, y′x)
            x, y = prex - y′x, prey + y′y
            return x ∈ 1..cols && y ∈ 1..rows ? interp(y, x) : z
        end
        next!(p)
    end
    rmat |> CTSinogram{M}
end


"""radon_alt(image::AbstractMatrix; <keyword arguments>)

Compute the Radon transform of `image` inside a circle of
radius `(min(rows,cols)-1) / 2` where `rows` and `cols` are the
dimensions of `image`.

See also: [`radon_std`](@ref)
"""
function radon_alt(
    image::AbstractMatrix{T};
    nd::Optional{Int} = nothing,
    nϕ::Optional{Int} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    ν::Real = 1,
    background::Optional{<:Real} = nothing,
    rescaled::Bool = true,
    interpolation::Optional{Interp} = nothing,
) where {T <: Real,Interp <: Union{Function,AbstractInterp2DOrNone}}
    M = typeof(image)
    rows, cols = size(image)
    l = min(rows, cols)
    nd = isnothing(nd) ? round(Int, l) : nd
    nϕ = isnothing(nϕ) ?
        2 * (round(Int, (rows + 1) * (cols + 1) / nd) ÷ 2) + 1 : nϕ
    x′₀ = (nd-l) ÷ 2
    t₀::T = (l - 1) / 2 * ν
    x₀::T = (cols - 1) / 2 + 1
    y₀::T = (rows - 1) / 2 + 1
    ϕ₀::T = deg2rad(α₀)
    Δϕ::T = deg2rad(α) / nϕ
    ts = linspace(-t₀..t₀, nd)
    scϕs = map(sincos, range(ϕ₀; step = Δϕ, length = nϕ))
    rimage = rescaled ? rescale(image) : image
    interp = isnothing(interpolation) ?
        interpolate(rimage) : interpolation(rimage)
    z::T = maybe(zero(T), background)
    rmat = fill(z, nd, nϕ)
    @info "Computing Radon transform..."
    p = Progress(nϕ, 0.2)
    Threads.@threads for iϕ ∈ 1:nϕ
        @inbounds s, c = scϕs[iϕ]
        @inbounds @simd for it ∈ 1:l
            t::T = ts[it+x′₀]
            prex::T = t * c + x₀
            prey::T = t * s + y₀
            for z ∈ ts
                x::T = prex - z * s
                y::T = prey + z * c
                if x ∈ 1..cols && y ∈ 1..rows
                    rmat[x′₀ + it,iϕ] += interp(y,x)
                end
            end
        end
        next!(p)
    end
    rmat |> CTSinogram{M}
end


const _default_radon_ref = Ref{Function}(radon_std)
default_radon() = _default_radon_ref[]
default_radon(other::Function) = _default_radon_ref[] = other
default_radon(other::Symbol) = @eval _default_radon_ref[] = $other


radon(image::AbstractMatrix; kwargs...) = default_radon()(image; kwargs...)


function radon(
    image::AbstractMatrix,
    geometry::AbstractParallelBeamGeometry;
    kwargs...
)
    nd = geometry.nd
    nϕ = geometry.nϕ
    α = geometry.α
    α₀ = geometry.α₀
    radon(image; nd, nϕ, α, α₀, kwargs...)
end


function radon(
    image::AbstractMatrix,
    geometry::AbstractFanBeamGeometry;
    kwargs...
)
    g′ = ParallelBeamGeometry(geometry)
    sinog = radon(image, g′; kwargs...)
    _, fan_sinog = para2fan(sinog, geometry)
    fan_sinog
end


radon(; kwargs...) = x -> radon(x; kwargs...)
radon(g::AbstractGeometry; kwargs...) = x -> radon(x, g; kwargs...)


struct Radon{Geometry <: AbstractGeometry} <: AbstractProjectionAlgorithm
    geometry::Geometry
end


datatype(rdn::Radon) = datatype(rdn.geometry)
alg_geometry(rdn::Radon) = rdn.geometry
alg_params(::Radon) = nothing


radon(image::AbstractMatrix, rdn::Radon; kwargs...) = radon(image, rdn.geometry; kwargs...)


(rdn::Radon)(image::AbstractMatrix; kwargs...) = radon(rdn, image; kwargs...)
(rdn::Radon)(; kwargs...) = radon(rdn; kwargs...)


project_image(image::CTImage, rdn::Radon; kwargs...) = image ↣ radon(rdn; kwargs...)
