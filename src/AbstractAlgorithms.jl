module AbstractAlgorithms

export AbstractCTAlgorithm, AbstractProjectionAlgorithm, AbstractReconstructionAlgorithm
export AbstractParams
export alg_geometry, alg_params
export project_image, reconstruct_image
export radon, iradon


using ..Monads


abstract type AbstractCTAlgorithm end
abstract type AbstractProjectionAlgorithm <: AbstractCTAlgorithm end
abstract type AbstractReconstructionAlgorithm <: AbstractCTAlgorithm end

abstract type AbstractParams end


function alg_geometry end
function alg_params end

function radon end
function iradon end
function project_image end
function reconstruct_image end


radon(alg::AbstractProjectionAlgorithm; kwargs...) = x -> radon(x, alg; kwargs...)
radon(x::Maybe, alg::AbstractProjectionAlgorithm; kwargs...) =
    x ↣ radon(alg; kwargs...)

for nm ∈ (:iradon, :reconstruct_image)
    @eval begin
        $nm(alg::AbstractReconstructionAlgorithm; kwargs...) =
            x -> $nm(x, alg; kwargs...)
        $nm(x::Maybe, alg::AbstractReconstructionAlgorithm; kwargs...) =
            x ↣ $nm(alg; kwargs...)
    end
end

end # module
