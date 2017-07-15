#! /usr/bin/env julia

# This files contains the abstract modelling of the lattice boltzmann method.
# It specifies abstract types, which are used in the function/ type definitions, includes after this files.

# The following three type specify a LBM modell
abstract type Collision end
abstract type Streaming end

# Models for different lattice stars
abstract type VelocitySet end
abstract type Flow end  # What kind of flow

abstract type Grid end
abstract type AbstractPatch <: Grid end # For collecting blocks
abstract type AbstractBlock <: Grid end # For containing populations

abstract type Geometry end
abstract type GeomProperty <: Geometry end
abstract type Solid <: Geometry end
abstract type Virtual <: Geometry end

# Boundary conditions are modeled onto a LBM type, since the modify the LBM.
abstract type Direction end
abstract type North <: Direction end
abstract type South <: Direction end
abstract type West <: Direction end
abstract type East <: Direction end

abstract type _1D <: VelocitySet end
abstract type _2D <: VelocitySet end
abstract type _3D <: VelocitySet end

abstract type SingleDistFlow <: Flow end
abstract type Compressible <: SingleDistFlow end

abstract type Incompressible <: SingleDistFlow end
abstract type FullPeriodicStreaming <: Streaming end
abstract type InnerStreaming <: Streaming end

abstract type BGK <: Collision end

export
    Collision,
    Streaming,
    Particle,
    Grid, Velocity_Set,
    Collision,
    Streaming,
    Grid, VelocitySet,
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

