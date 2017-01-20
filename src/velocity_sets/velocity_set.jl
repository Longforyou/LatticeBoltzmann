#! /usr/bin/env julia

module velocity_set

using Abstract_LBM
"""
This module contains definitions of types for referencing
the different kinds of velocity fields.
Example
-------
using velocity_set._D2Q9

make the weights w and the direction arrays c_x & c_y
 visible in the global namespace
"""

abstract _2D <: Velocity_Set
abstract _3D <: Velocity_Set


# Include the 2D velocity sets
include("velocity_2d.jl"); export _D2Q9

# TODO 3d velocity sets

include("equilibrium_func.jl")
include("macro_var.jl")

end # module velocity_set
