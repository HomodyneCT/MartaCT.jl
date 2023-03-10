module FDKAlgorithm

Base.Experimental.@optlevel 3

export AbstractFDK
export FDK, FDKGPU
export FDKInfo
export fdk


using IntervalSets
using Unitful: NoUnits
using ..Interpolation: blerp
using ..Utils: linspace
using ..Geometry
import ..AbstractAlgorithms: reconstruct_image
using ..AbstractAlgorithms


abstract type AbstractFDK <: AbstractIRadonAlgorithm end
abstract type AbstractFDKData end


function fdkdata end
function fdkfilter end


include("fdk/fdk.jl")


struct FDKInfo{
    G <: AbstractConeBeamGeometry,
    A <: AbstractFDK
} <: AlgorithmInfo{A}
    geometry::G
    algorithm::A
end

FDKInfo(g::AbstractConeBeamGeometry) = FDKInfo(g, FDK())


using Requires

function __init__()
    @static if Sys.isapple()
        @debug "Requiring Metal"
        @require Metal="dde4c033-4e86-420c-a63e-0dd931031962" include("fdk/fdk_gpu.jl")
    end
end

end # module
