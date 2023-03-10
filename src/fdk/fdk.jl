import ..Filters: AbstractCTFilter, RamLak, CTFilterPlan

struct FDK{F<:AbstractCTFilter} <: AbstractFDK end

FDK(::Type{F} = RamLak) where {F<:AbstractCTFilter} = FDK{F}()


function fdk end


struct FDKData{P<:AbstractArray,P′<:AbstractArray,W<:AbstractMatrix,A<:AbstractArray,FP<:CTFilterPlan} <: AbstractFDKData
    sz::NTuple{3,Int}
    proj::P
    fproj::P′
    weights::W
    vdata::A
    filter_plan::FP
end


function fdkdata(::Type{F}, proj, us, vs, xs, ys, βs) where {F<:AbstractCTFilter}
    T = promote_type(eltype(us), eltype(vs))
    nv, nu, nβ = size(proj)
    @assert length(us) == nu && length(vs) == nv
    M = zeros(eltype(proj), nextprod((2,), nu), nv, nβ)
    fp = CTFilterPlan(F, M, 1)
    P = plan(fp)
    M̂ = zeros(complex(eltype(M)), axes(inv(P)))
    _1 = oneunit(T)
    W = @. inv(sqrt(_1 + us^2 + vs'^2))
    vdata = _vdata(xs′, ys′, βs)
    FDKData((nu,nv,nβ), M, M̂, W, vdata, fp)
end


function fdkfilter(proj, fdk_data)
    nu = fdk_data.sz[1]
    M = fdk_data.proj
    M̂ = fdk_data.fproj
    W = fdk_data.weights
    fp = fdk_data.filter_plan
    M′ = @view M[begin:begin+nu-1,:,:]
    permutedims!(M′, proj, (2,1,3))
    @. M′ = W * M′
    M[begin+nu:end,:,:] .= zero(eltype(M))
    reverse!(M, dims=2)
    fill!(M̂, zero(eltype(M̂)))
    ctfilter!(M̂, M, fp)
    M′
end


@inline function _compute_backprojection!(tomog, M, matV, xs, ys, zs, βs, u₀, v₀, δu, δv)
    T = eltype(tomog)
    _0 = zero(T)
    γ::T = inv(4π^2)
    nu, nv, _ = size(M)
    u₀ += oneunit(u₀)
    v₀ += oneunit(v₀)
    iδu, iδv = inv(δu), inv(δv)
    xs′ = @. T(NoUnits(xs * iδu))
    ys′ = @. T(NoUnits(ys * iδu))
    zs′ = @. T(NoUnits(zs * iδv))
    Threads.@threads for ind in CartesianIndices(tomog)
        j, i, k = ind.I
        @inbounds begin
            x, y, z = xs′[i], ys′[j], zs′[k]
            t = _0
            @simd for b in eachindex(βs)
                σ = matV[j, i, b]
                s, c = βs[b]
                u′ = u(x, y, s, c, σ) + u₀
                v′ = v(z, σ) + v₀
                if u′ ∈ 1..nu && v′ ∈ 1..nv
                    t += σ^2 * blerp(M, u′, v′, b)
                end
            end
            tomog[ind] = t * γ
        end
    end
    tomog
end


function fdk(::Type{F}, proj, xs::AbstractVector, ys::AbstractVector, zs::AbstractVector; R, px=nothing, δ=nothing) where {F<:AbstractCTFilter}
    T = eltype(proj)
    δu, δv = something(δ, px, (oneunit ∘ first).((xs, zs))) |> _make_dpx
    nv, nu, nβ = size(proj)
    u₀::T, v₀::T = @. ((nu, nv) - 1) / 2
    us = range(-u₀..u₀, nu) * T(δu / R)
    vs = range(-v₀..v₀, nv) * T(δv / R)
    βs = sincospi.(linspace(T, 0..2(nβ-1)//nβ, nβ))
    @debug "Filtering"
    fdk_data = fdkdata(F, proj, us, vs, xs / R, ys / R, βs)
    M = fdkfilter(proj, fdk_data)
    @debug "V matrix"
    VV = fdk_data.vdata
    @debug "Backprojection"
    tomog = similar(VV, length.((xs, ys, zs)))
    _compute_backprojection!(tomog, M, VV, xs, ys, zs, βs, u₀, v₀, δu, δv)
end


function _vdata(::Type{T}, xs, ys, βs) where T
    @assert unit(eltype(xs)) == unit(eltype(ys))
    nx, ny, nb = length.((xs, ys, βs))
    vdata = Array{T,3}(undef, ny, nx, nb)
    xs′ = @. T(NoUnits(xs))
    ys′ = @. T(NoUnits(ys))
    βs′ = @. T(βs)
    _vdata!(vdata, xs′, ys′, βs′)
end


function _vdata!(vdata, xs, ys, βs)
    Threads.@threads for ind in CartesianIndices(vdata)
        j, i, k = ind.I
        @inbounds x, y, (s, c) = xs′[i], ys′[j], βs[k]
        @inbounds vdata[ind] = V(x, y, s, c)
    end
    vdata
end
