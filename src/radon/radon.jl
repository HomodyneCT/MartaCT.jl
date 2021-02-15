import LoopVectorization: @avx, vmapreduce

function radon(
    image::AbstractMatrix{T};
    nd::Optional{Int} = nothing,
    nϕ::Optional{Int} = nothing,
    α::Real = 360,
    α₀::Real = 0,
    interpolation::Optional{Interp} = nothing,
) where {T <: Real, Interp <: Union{Function,AbstractBilinearInterpolation}}
    M = typeof(image)
    rows, cols = size(image)
    nd = maybe(round(Int, sqrt(rows * rows + cols * cols)), nd)
    nϕ = maybe(2 * (round(Int, (rows + 1) * (cols + 1) / nd) ÷ 2) + 1, nϕ)

    x′₀::T = T(nd - 1) / 2
    x₀::T = T(cols + 1) / 2
    y₀::T = T(rows + 1) / 2
    t₀ = (x₀, y₀)
    ϕ₀::T = deg2rad(α₀)
    Δϕ::T = deg2rad(α) / nϕ

    x′s = range(-1, 1; length = nd) |> collect
    # x′s = range(-1; step = 2 / nd, length = nd) |> collect
    # ϕs = map(reim ∘ cis, (0:nϕ - 1) * Δϕ .+ ϕ₀)
    ϕs = map(reim ∘ cis, range(ϕ₀; step = Δϕ, length = nϕ))

    rmat = zeros(T, nd, nϕ)
    indices = Vector{NTuple{2,Int}}(undef, length(rmat))
    x′ϕs = Vector{NTuple{2,T}}(undef, length(rmat))

    Threads.@threads for k in LinearIndices(rmat)
        ix = (k - 1) % nd
        iϕ = (k - 1) ÷ nd
        @inbounds indices[k] = k, iϕ
        @inbounds x′ϕs[k] = x′s[ix + 1] .* ϕs[iϕ + 1]
    end

    rimage = rescale(image)
    interp = isnothing(interpolation) ? interpolate(rimage) : interpolation(rimage)

    # interp_value = (y, x) -> begin
    #     1 ≤ x ≤ cols && 1 ≤ y ≤ rows && return interp(y, x)
    #     return 0
    # end

    compute_projection = (k, iϕ) -> begin
        @inbounds x′x, x′y = x′ϕs[k]
        prex = (x′x + 1) * x₀
        prey = (x′y + 1) * y₀
        offset = iϕ * nd
        #@views sum(x′ϕs[offset + 1:offset + nd]) do (y′y, y′x)
        #s::T = zero(T)
        # @avx for ii ∈ view(x′ϕs, offset + 1:offset + nd)
        #     (y′y, y′x) = ii
        #     x = prex - y′x * x₀
        #     y = prey + y′y * y₀
        #     if 1 ≤ x ≤ cols && 1 ≤ y ≤ rows
        #         s += interp(y, x)
        #     end
        # end
        # return s
        vmapreduce(+, view(x′ϕs, offset + 1:offset + nd)) do (y′y, y′x)
            x = prex - y′x * x₀
            y = prey + y′y * y₀
            if 1 ≤ x ≤ cols && 1 ≤ y ≤ rows
                return interp(y, x)
            end
        end
    end

    @info "Computing Radon transform..."
    p = Progress(length(rmat), 0.2)

    Threads.@threads for (k, iϕ) in indices
        @inbounds rmat[k] = compute_projection(k, iϕ)
        next!(p)
    end

    rmat |> CTSinogram{M}
end


# function radon(
#     image::AFMatrix{T};
#     nd = nothing,
#     nϕ = nothing,
#     α = 360,
#     α₀ = 0,
#     interpolation = nothing,
# ) where {T<:Real}
#     rows, cols = size(image)
#     if isnothing(nd)
#         nd = round(Int, sqrt(rows * rows + cols * cols))
#     end
#     if isnothing(nϕ)
#         nϕ = 2 * (round(Int, (rows + 1) * (cols + 1) / nd) ÷ 2) + 1
#     end
#     if isnothing(interpolation)
#         interpolation = mat -> LinearInterpolation(axes(mat), mat)
#     end
#     x′₀::T = T(nd - 1) / 2
#     x₀::T = T(cols + 1) / 2
#     y₀::T = T(rows + 1) / 2
#     y′s = (0:nd-1) .- x′₀
#     Δϕ::T = deg2rad(α) / nϕ
#     ϕ₀::T = deg2rad(α₀)
#     rmat = AFArray(zeros(T, nd, nϕ))
#     m = length(rmat)
#     interp = interpolation(image)
#
#     @info Computing Radon transform..."
#
#     for k = 1:m
#         ϕ::T = (((k - 1) ÷ nd) * Δϕ + ϕ₀)
#         x′ = (k - 1) % nd - x′₀
#         cϕ::T = cos(ϕ)
#         sϕ::T = sin(ϕ)
#         prex = x′ * cϕ + x₀
#         prey = x′ * sϕ + y₀
#         sum = zero(T)
#         for y′ in y′s
#             x = prex - y′ * sϕ
#             y = prey + y′ * cϕ
#             if 1 <= x <= cols && 1 <= y <= rows
#                 sum += interp(y, x)
#             end
#         end
#         rmat[k] = sum
#     end
#
#     rmat
# end


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
