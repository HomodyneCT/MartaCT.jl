module AbstractAlgorithms

export
    AbstractCTAlgorithm,
    AbstractProjectionAlgorithm,
    AbstractReconstructionAlgorithm,
    AbstractIRadonAlgorithm
export AbstractParams
export AlgorithmInfo
export project_image, reconstruct_image
export radon, iradon
export alg_name, alg_params


using SimpleTraits
using ..Monads
using ..Geometry:
    AbstractGeometry, AbstractParallelBeamGeometry, AbstractFanBeamGeometry
using ..FanBeam: fan2para
using ..Coordinates


@inline function _alg_progress(
    ::Type{T},
    desc::AbstractString,
    n::Integer,
    enabled::Bool,
    dt::Real=0.2
) where T
    T(n; dt, desc, enabled)
end


abstract type AbstractCTAlgorithm end
abstract type AbstractProjectionAlgorithm <: AbstractCTAlgorithm end
abstract type AbstractReconstructionAlgorithm <: AbstractCTAlgorithm end
abstract type AbstractIRadonAlgorithm <: AbstractReconstructionAlgorithm end
abstract type AlgorithmInfo{A<:AbstractCTAlgorithm} end
abstract type AbstractParams end

alg_name(::Type{A}) where {A <: AbstractCTAlgorithm} = nameof(A)
alg_name(a::AbstractCTAlgorithm) = alg_name(typeof(a))
alg_name(ainfo::AlgorithmInfo) = alg_name(ainfo.algorithm)
@inline function alg_params(ainfo::AInfo) where {AInfo <: AlgorithmInfo}
    :params ∈ fieldnames(AInfo) && return Some(ainfo.params)
    return nothing
end


function radon end
function iradon end
function project_image end
function reconstruct_image end


const _alg_fns =
    (:radon, :iradon, :project_image, :reconstruct_image)
const _algfn_map = (;
    :AbstractProjectionAlgorithm => (:radon, :project_image),
    :AbstractReconstructionAlgorithm => (:iradon, :reconstruct_image),
)

for (T, fns) ∈ pairs(_algfn_map)
    for nm ∈ fns
        @eval begin
            @inline $nm(; kwargs...) = x -> $nm(x; kwargs...)
            @inline $nm(alg::$T; kwargs...) = x -> $nm(x, alg; kwargs...)
            @traitfn @inline function $nm(
                x::M, alg::$T; kwargs...
            ) where {M; Monad{M}}
                x ↣ $nm(alg; kwargs...)
            end
            @traitfn @inline function $nm(x::M; kwargs...) where {M; Monad{M}}
                x ↣ $nm(; kwargs...)
            end
            @inline function $nm(
                ainfo::AlgorithmInfo{A}; kwargs...
            ) where {A<:$T}
                x -> $nm(x, ainfo; kwargs...)
            end
            @inline function $nm(
                x::AbstractMatrix,
                ainfo::AlgorithmInfo{A};
                kwargs...
            ) where {A<:$T}
                g = ainfo.geometry
                p = alg_params(ainfo)
                a = ainfo.algorithm
                isnothing(p) && return $nm(x, g, a; kwargs...)
                $nm(x, g, something(p), a; kwargs...)
            end
            @traitfn @inline function $nm(
                x::M,
                ainfo::AlgorithmInfo{A};
                kwargs...
            ) where {M,A<:$T; Monad{M}}
                x ↣ $nm(ainfo; kwargs...)
            end
        end
    end
end


macro _defapi(f::Symbol, g::Symbol, A::Symbol)
    esc(quote
        @inline function $f(
            image::AbstractMatrix,
            alg::$A;
            kwargs...
        )
            $g(image, alg; kwargs...)
        end
        @inline function $f(
            image::AbstractMatrix,
            geometry::AbstractGeometry,
            alg::$A;
            kwargs...
        )
            $g(image, geometry, alg; kwargs...)
        end
        @inline function $f(
            image::AbstractMatrix,
            geometry::AbstractGeometry,
            params::AbstractParams,
            alg::$A;
            kwargs...
        )
            $g(image, geometry, params, alg; kwargs...)
        end
    end)
end


@_defapi project_image radon AbstractProjectionAlgorithm
@_defapi reconstruct_image iradon AbstractReconstructionAlgorithm


@inline function iradon(
    sinog::AbstractMatrix,
    alg::AbstractReconstructionAlgorithm;
    kwargs...
)
    alg(sinog; kwargs...)
end

@inline function iradon(
    sinog::AbstractMatrix,
    xs::AbstractVector,
    ys::AbstractVector,
    alg::AbstractReconstructionAlgorithm,
    coo::AbstractCoordinates = Cartesian();
    kwargs...
)
    alg(sinog, xs, ys, coo; kwargs...)
end

@inline function iradon(
    sinog::AbstractMatrix,
    xs::AbstractVector,
    ys::AbstractVector,
    alg::AbstractReconstructionAlgorithm,
    coo::AbstractCoordinates = Cartesian();
    kwargs...
)
    alg(sinog, xs, ys, ϕs, coo; kwargs...)
end

@inline function iradon(
    sinog::AbstractMatrix,
    xs::AbstractVector,
    alg::AbstractReconstructionAlgorithm,
    coo::AbstractCoordinates = Cartesian();
    kwargs...
)
    alg(sinog, xs, coo; kwargs...)
end

@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractParallelBeamGeometry,
    alg::AbstractReconstructionAlgorithm;
    kwargs...
)
    rows = geometry.rows
    cols = geometry.cols
    α = geometry.α
    α₀ = geometry.α₀
    iradon(sinog, alg; rows, cols, α, α₀, kwargs...)
end

@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractFanBeamGeometry,
    alg::AbstractReconstructionAlgorithm;
    kwargs...
)
    g′, para_sinog = fan2para(sinog, geometry)
    iradon(para_sinog, g′, alg; kwargs...)
end

@inline function iradon(
    g::AbstractGeometry,
    alg::AbstractReconstructionAlgorithm;
    kwargs...
)
    x -> iradon(x, g, alg; kwargs...)
end

@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractParallelBeamGeometry,
    params::AbstractParams,
    alg::AbstractReconstructionAlgorithm;
    kwargs...
)
    rows = geometry.rows
    cols = geometry.cols
    α = geometry.α
    α₀ = geometry.α₀
    alg(sinog, params; rows, cols, α, α₀, kwargs...)
end


@inline function iradon(
    sinog::AbstractMatrix,
    geometry::AbstractFanBeamGeometry,
    params::AbstractParams,
    alg::AbstractReconstructionAlgorithm;
    kwargs...
)
    g′, para_sinog = fan2para(sinog, geometry)
    iradon(sinog, g′, params, alg; kwargs...)
end

end # module
