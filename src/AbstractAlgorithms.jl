module AbstractAlgorithms

export
    AbstractCTAlgorithm,
    AbstractProjectionAlgorithm,
    AbstractReconstructionAlgorithm,
    AbstractIRadonAlgorithm
export AbstractParams
export project_image, reconstruct_image
export radon, iradon


using SimpleTraits
using ..Monads


abstract type AbstractCTAlgorithm end
abstract type AbstractProjectionAlgorithm <: AbstractCTAlgorithm end
abstract type AbstractReconstructionAlgorithm <: AbstractCTAlgorithm end
abstract type AbstractIRadonAlgorithm <: AbstractReconstructionAlgorithm end

abstract type AbstractParams end


function radon end
function iradon end
function project_image end
function reconstruct_image end


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
        end
    end
end


macro _alg_progress(desc, n, p, dt=0.2, T=:Progress)
    :($T($n; dt=$dt, desc=$desc, enabled=$p))
end

end # module
