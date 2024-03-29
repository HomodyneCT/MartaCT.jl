img = SquareImage(rows=4, cols=4, l=2, calibration_value=1, background=0)
@testset "Testing FBP algorithms" begin
    @testset "Testing FBP FFT algorithm" begin
        sinog = radon(img, Radon())
        tomog = iradon(sinog, FBP())
        @test tomog ≈ [
            [-0.642772  0.374742  0.374742  -0.642772];
            [0.033042  1.05056   1.05056    0.0330419];
            [0.033042  1.05056   1.05056    0.0330422];
            [-0.642772  0.374742  0.374742  -0.642772]
        ]
    end
    @testset "Testing FBP FFT square algorithm" begin
        sinog = radon(img, RadonSquare())
        tomog = iradon(sinog, FBPFFTSquare())
        @test tomog ≈ [
            [0.364821   0.55048   0.55428  -0.0358717];
            [0.270617   1.22652   1.21234   0.0387748];
            [0.0387749  1.21234   1.22652   0.270617];
            [-0.0358716  0.554279  0.55048   0.364821]
        ]
    end
end
