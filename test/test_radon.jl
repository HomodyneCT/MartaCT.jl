img = SquareImage(rows=4, cols=4, l=2, calibration_value=1, background=0)
@testset "Testing Radon transform" begin
    @test img ≈ [
        [0 0 0 0];
        [0 1 1 0];
        [0 1 1 0];
        [0 0 0 0]
    ]
    @testset "Testing diag algorithm" begin
        sinog = radon(img, Radon())
        @test sinog ≈ [
            [0.0       0.0       0.0       0.0       0.0];
            [0.178452  0.145332  0.164798  0.164798  0.145332];
            [0.785413  0.763186  0.745526  0.745525  0.763186];
            [0.785413  0.763186  0.745526  0.745525  0.763186];
            [0.178452  0.145332  0.164798  0.164798  0.145332];
            [0.0       0.0       0.0       0.0       0.0]
        ]
    end
    @testset "Testing square algorithm" begin
        sinog = radon(img, RadonSquare())
        @test sinog ≈ [
            [0.0       0.05596    0.044909  0.064763   0.0647629  0.044909  0.05596];
            [0.888889  0.764706   0.846043  0.804828   0.804828   0.846043  0.764706];
            [0.888889  0.764706   0.846043  0.804828   0.804828   0.846043  0.764706];
            [0.0       0.0559601  0.044909  0.0647629  0.0647629  0.044909  0.05596]
        ]
    end
end
