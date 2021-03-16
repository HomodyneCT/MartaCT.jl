module CalibrationBase

Base.Experimental.@optlevel 3

using SimpleTraits
using ..Monads


function calibration_data end
function calibrate_image end
function calibrate_tomogram end


for nm ∈ (:calibrate_image, :calibrate_tomogram)
    @eval begin
        $nm(; kwargs...) = x -> $nm(x; kwargs...)
        @traitfn $nm(m::M; kwargs) where {M; Monad{M}} =
            mmap($nm(; kwargs...), m)
        @traitfn $nm(m::M, x; kwargs) where {M; Monad{M}} =
            mmap($nm(x; kwargs...), m)
    end
end

end # module
