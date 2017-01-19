#! /usr/bin/env julia

module velocity_set

# definitions of types for referencing constant field etc.
abstract Veloctiy_Set

export Veloctiy_Set

# Include the 2D velocity sets
include("velocity_2d.jl"); export _D2Q9

# TODO 3d velocity sets

end # module velocity_set
