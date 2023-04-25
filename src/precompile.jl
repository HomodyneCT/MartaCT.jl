using PrecompileTools

@setup_workload begin
    _pre_types = (:Float32, :Float64)
    width = 5
    nϕ = 11

    @compile_workload begin
        for T ∈ _pre_types, Img ∈ CTTestImages._gs_images, G ∈ Geometry._geometry_names, A ∈ (:Radon, :RadonSquare), F ∈ (:FBP, :FBPAFFT, :FBPAFFTSquare, :FBPFFTSquare)
            @eval begin
                gs = $Img($T; width=$width)
                g = $G(gs; nϕ=$nϕ)
                ct = project_image(gs, g, $A()) |> reconstruct_image($F())
                ctc = calibrate_tomogram(ct, gs)
            end
        end
    end
end
