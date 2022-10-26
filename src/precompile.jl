using SnoopPrecompile

@precompile_setup begin
    const _pre_types = (:Float32, :Float64)
    const width = 5
    const nϕ = 5

    @precompile_all_calls begin
        for T ∈ _pre_types, Img ∈ CTTestImages._gs_images, G ∈ Geometry._geometry_names
            @info "Executing... " T Img G
            @eval begin
                gs = $Img($T; width)
                g = $G(gs; nϕ)
                ct = FBPScanner(g, gs) |> project_image |> reconstruct_image
                ctc = calibrate_tomogram(ct, gs)
            end
        end
    end
end
