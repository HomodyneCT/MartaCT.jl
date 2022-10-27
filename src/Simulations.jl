module Simulations

Base.Experimental.@optlevel 3

export AbstractSimulation, CTSimulation
export simulate

using ..Monads
using ..Utils
using ProgressMeter, LinearAlgebra, IntervalSets
import Random, Distributions, StatsBase

abstract type AbstractSimulation end

function simulate end
function simulate_blocks end
function simulate_single end

include("simulations/sinogram_sampling.jl")


@inline function Random.rand(
    rng::Random.AbstractRNG,
    s::Random.SamplerTrivial{ORI{T}},
) where {T <: AbstractFloat}
    a, b = endpoints(s[])
    t = rand(rng, T)
    a * (one(T) - t) + b * t
end

end # module
