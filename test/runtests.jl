include("setup.jl")

tests = ["radon", "fbp", "geometry"]

if !isempty(ARGS)
    tests = ARGS
end

@testset "Tests for the MartaCT package" begin
    for t ∈ tests
        include("test_$t.jl")
    end
end
