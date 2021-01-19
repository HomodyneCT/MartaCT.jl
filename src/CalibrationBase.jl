module CalibrationBase

Base.Experimental.@optlevel 3

using ..Monads


function calibrate_image end
function calibrate_tomogram end


for nm ∈ (:calibrate_image, :calibrate_tomogram)
    @eval begin
        $nm(; kwargs...) = x -> $nm(x; kwargs...)
        $nm(m::Maybe; kwargs) = m ↣ $nm(; kwargs...)
    end
end

end # module
