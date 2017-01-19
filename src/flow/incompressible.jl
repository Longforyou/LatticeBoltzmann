#! /usr/bin/env julia

abstract Incompressible <: Flow

# Type definition for the lattice boltzmann 2D mesh
type LBM_Incompressible{V, F, S, C} <: LBM{V<:Velocity_Set,
                                           F<:Incompressible,
                                           S<:Streming,
                                           C<:Collision}
  
  constants::LBM_Constants
  grid::Grid{T}
  bound::Array{Boundary, 1} # For Bounces etc...

  Lattice_Boltzmann_2D(constants::LBM_Constants,
                       grid::Grid{T},
                       bound::Array{Boundary, 1}) =
    new(constants, grid, bound)

end
