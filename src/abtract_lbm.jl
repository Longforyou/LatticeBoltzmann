#! /usr/bin/env julia

"""
This files contains the abstract modelling of the lattice boltzmann method.
It specifies abstract types, which are used in the function/ type definitions, includes after this files.
"""

# The following three type specify a LBM modell
abstract Collision 
abstract Streaming
abstract Velocity_Set
abstract Flow

abstract LBM{V<:Velocity_Set, F <: Flow, S <: Streaming, C <: Collision}

# Boundary conditions are modeled onto a LBM type, since the modify the LBM.
