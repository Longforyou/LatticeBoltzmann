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

abstract _1D <: Velocity_Set
abstract _2D <: Velocity_Set
abstract _3D <: Velocity_Set

abstract Compressible <: Flow
abstract Incompressible <: Flow
abstract FullPeriodicStreaming <: Streaming
abstract InnerStreaming <: Streaming

abstract BGK <: Collision

export
    LBM,
    Collision,
    Streaming,
    Velocity_Set,
    _1D,_2D, _3D,
    Flow, Direction,
    North, South,
    West, East,
    Compressible, Incompressible,
    FullPeriodicStreaming,
    InnerStreaming, BGK

end # module Abstract_LBM
