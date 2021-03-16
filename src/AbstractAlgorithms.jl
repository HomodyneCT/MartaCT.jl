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
            $nm(; kwargs...) = x -> $nm(x; kwargs...)
            $nm(alg::$T; kwargs...) = x -> $nm(x, alg; kwargs...)
            @traitfn function $nm(x::M, alg::$T; kwargs...) where {M; Monad{M}}
                x ↣ $nm(alg; kwargs...)
            end
            @traitfn function $nm(x::M; kwargs...) where {M; Monad{M}}
                x ↣ $nm(; kwargs...)
            end
            function $nm(ainfo::AlgorithmInfo{A}; kwargs...) where {A<:$T}
                x -> $nm(x, ainfo; kwargs...)
            end
            function $nm(
                x::AbstractMatrix,
                ainfo::AlgorithmInfo{A};
                kwargs...
            ) where {A<:$T}
                g = ainfo.geometry
                p = alg_params(ainfo)
                a = ainfo.algorithm
                isnothing(p) && return a(x, g; kwargs...)
                a(x, g, something(p); kwargs...)
            end
            @traitfn function $nm(
                x::M,
                ainfo::AlgorithmInfo{A};
                kwargs...
            ) where {M,A<:$T; Monad{M}}
                x ↣ $nm(ainfo; kwargs...)
            end
        end
    end
end


macro _alg_progress(desc, n, p, dt=0.2, T=:Progress)
    :($T($n; dt=$dt, desc=$desc, enabled=$p))
end

end # module
