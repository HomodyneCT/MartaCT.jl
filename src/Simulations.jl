module Simulations

using ..Monads
using ..CTImages
using ..Utils: linspace
using ProgressMeter, LinearAlgebra
import Random, Distributions, StatsBase
const Rnd = Random
const Dist = Distributions
const SB = StatsBase

abstract type AbstractSimulation end

include("simulations/sinogram_sampling.jl")

end # module
