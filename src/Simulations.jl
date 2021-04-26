module Simulations

Base.Experimental.@optlevel 3

export AbstractSimulation, CTSimulation
export simulate

using ..Monads
using ..CTImages
using ..Utils: linspace, half
using ProgressMeter, LinearAlgebra, IntervalSets
import Random, Distributions, StatsBase

abstract type AbstractSimulation end

function simulate end

include("simulations/sinogram_sampling.jl")

end # module
