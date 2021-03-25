module Coordinates

export AbstractCoordinates
export Cartesian, Polar

abstract type AbstractCoordinates end

struct Cartesian <: AbstractCoordinates end
struct Polar <: AbstractCoordinates end

end # module
