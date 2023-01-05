module FDKAlgorithm

Base.Experimental.@optlevel 3

export AbstractFDK
export FDK, FDKGPU
export FDKInfo

using FFTW, LinearAlgebra, IntervalSets
using Unitful: NoUnits
using ..Interpolation: blerp
using ..Utils: linspace
import ..AbstractAlgorithms: reconstruct_image
using ..AbstractAlgorithms, ..Coordinates


include("fdk/fdk.jl")


using Requires

function __ini__()
    @static if Sys.isapple()
        @debug "Requiring Metal"
        @require Metal="dde4c033-4e86-420c-a63e-0dd931031962" include("fdk/fdk_gpu.jl")
    end
end

end # module
