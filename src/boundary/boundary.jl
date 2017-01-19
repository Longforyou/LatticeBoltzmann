#! /usr/bin/env julia

#= 
Basic idea is provide indices during the computations.
This provides other function to directly compute the data

=#

# Abstract type for the dispatch
abstract SolType
abstract Eq <: SolType
abstract NonEq <: SolType
abstract NonEqBounce <: SolType
abstract BounceCondition <: Boundary

include("neumann.jl")
include("dirichlet.jl")
include("open.jl")
include("bounce.jl")
include("corner.jl")
include("periodic_dirichlet.jl")
