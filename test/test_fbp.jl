img = SquareImage(rows=4, cols=4, l=2, calibration_value=1, background=0)
@testset "Testing FBP algorithms" begin
    @testset "Testing FBP FFT algorithm" begin
        sinog = radon(img, Radon())
        tomog = iradon(sinog, FBP())
        @test tomog ≈ [
            [0.693953  1.74108    1.48404    0.357039];
            [0.708007  1.75514    1.49809    0.371093];
            [-0.24023   0.806901   0.549855  -0.577144];
            [-0.926158  0.120973  -0.136073  -1.26307]
        ]
    end
    @testset "Testing FBP FFT square algorithm" begin
        sinog = radon(img, RadonSquare())
        tomog = iradon(sinog, FBPFFTSquare())
        @test tomog ≈ [
            [0.177991   0.55048  0.554279  -0.260271];
            [0.0837881  1.22652  1.21234   -0.185625];
            [-0.148054   1.21234  1.22652   -0.132705];
            [-0.0358718  0.55428  0.55048    0.140421]
        ]
    end
end
