#! /usr/bin/env julia

# This files contains the abstract modelling of the lattice boltzmann method.
# It specifies abstract types, which are used in the function/ type definitions, includes after this files.

# The following three type specify a LBM modell
abstract Collision 
abstract Streaming
<<<<<<< HEAD
abstract Velocity_Set # Models for different lattice stars
=======
abstract VelocitySet # Models for different lattice stars
abstract PopulationSet
>>>>>>> origin/dev2
abstract Flow  # What kind of flow 
abstract AbstractDomain
abstract AbstractBlock
abstract Geometry
abstract GeomProperty <: Geometry
abstract Solid <: Geometry
abstract Virtual <: Geometry

# Boundary conditions are modeled onto a LBM type, since the modify the LBM.
abstract Direction
abstract North <: Direction
abstract South <: Direction
abstract West <: Direction
abstract East <: Direction


# Grid is seperated into the two parts. The individual lattice and the grid which holds all lattices
abstract Particle

abstract Grid
abstract _1D <: VelocitySet
abstract _2D <: VelocitySet
abstract _3D <: VelocitySet

abstract SingleDistFlow <: Flow
abstract Compressible <: SingleDistFlow
abstract Incompressible <: SingleDistFlow
abstract FullPeriodicStreaming <: Streaming
abstract InnerStreaming <: Streaming

<<<<<<< HEAD
export
    Collision,
    Streaming,
    Particle,
    Grid, Velocity_Set,
=======
abstract SinglePopulationSet <: PopulationSet

abstract BGK <: Collision

export
    Collision,
    Streaming,
    Grid, VelocitySet,
>>>>>>> origin/dev2
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

