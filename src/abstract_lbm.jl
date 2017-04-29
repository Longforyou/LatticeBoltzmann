#! /usr/bin/env julia

"""
This files contains the abstract modelling of the lattice boltzmann method.
It specifies abstract types, which are used in the function/ type definitions, includes after this files.
"""
module Abstract_LBM

# The following three type specify a LBM modell
abstract Collision 
abstract Streaming
abstract Velocity_Set # Models for different lattice stars
abstract Population_Set
abstract Flow  # What kind of flow 
abstract AbstractDomain
abstract AbstractBlock

# Boundary conditions are modeled onto a LBM type, since the modify the LBM.
abstract Direction
abstract North <: Direction
abstract South <: Direction
abstract West <: Direction
abstract East <: Direction


# Grid is seperated into the two parts. The individual lattice and the grid which holds all lattices
abstract Particle

abstract Grid
abstract _1D <: Velocity_Set
abstract _2D <: Velocity_Set
abstract _3D <: Velocity_Set

abstract Compressible <: Flow
abstract Incompressible <: Flow
abstract FullPeriodicStreaming <: Streaming
abstract InnerStreaming <: Streaming

abstract BGK <: Collision

export
    Collision,
    Streaming,
    Particle,
    Grid, Velocity_Set,
    AbstractBlock,
    AbstractDomain,
    Population_Set,
    _1D,_2D, _3D,
    Flow, Direction,
    North, South,
    West, East,
    Compressible, Incompressible,
    FullPeriodicStreaming,
    InnerStreaming, BGK

end # module Abstract_LBM
