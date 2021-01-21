using Marta

const _pre_types = (:Float32, :Float64)
const width = 5
const nϕ = 5

for T ∈ _pre_types, Img ∈ TestImages._gs_images, G ∈ Geometry._geometry_names
    @eval begin
        gs = $Img($T; width)
        FBPScanner(gs)
        g = $G(gs; nϕ)
        ct = FBPScanner(g, gs) |> project_image |> reconstruct_image
        ctc = calibrate_tomogram(ct, gs)
    end
end
