using Marta

const _pre_types = (:Float32, :Float64)
const _pre_imgs = (:GrayScaleLine, :GrayScalePyramid, :WhiteRect)
const _pre_geoms = (:ParallelBeamGeometry, :FanBeamGeometry)
const width = 5
const nϕ = 5

for t ∈ _pre_types, img ∈ _pre_imgs, g ∈ _pre_geoms
    @eval begin
        T = $t
        gs = $img(T; width)
        pbg = $g(T; nϕ, gs.width, gs.height)
        ct = FBPScanner(pbg, gs) |> project_image |> reconstruct_image
        ctc = calibrate_tomogram(ct, gs)
    end
end
