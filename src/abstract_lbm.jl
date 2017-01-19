#! /usr/bin/env julia

"""
This files contains the abstract modelling of the lattice boltzmann method.
It specifies abstract types, which are used in the function/ type definitions, includes after this files.
"""
module Abstract_LBM

# The following three type specify a LBM modell
abstract Collision 
abstract Streaming
abstract Velocity_Set
abstract Flow

abstract LBM{V <: Velocity_Set, F <: Flow,
             S <: Streaming, C <: Collision}

# Boundary conditions are modeled onto a LBM type, since the modify the LBM.
abstract Direction
abstract North <: Direction
abstract South <: Direction
abstract West <: Direction
abstract East <: Direction

export
    Collision,
    Streaming,
    Velocity_Set,
    Flow, Direction,
    North, South,
    West, East

end # module Abstract_LBM
