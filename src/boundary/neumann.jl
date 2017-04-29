#! /usr/bin/env julia

using ._D2Q9

# ====== Neumann conditions

"""
    Neumann{T<:Direction, S<:SolType, V<:Velocity_Set}

Definition of a Neumann-Boundary-Condition. Contains three
parametric types for dispatching on different implementations.
TODO implement differeny SolTypes
"""
immutable Neumann{T<:Direction, S<:SolType,
                  V <: Velocity_Set} <: Boundary
  vel::Float64
  rows::Array{Int64, 1}
  cols::Array{Int64, 1}

end

"""
    compute_neumann{T<:Direction}(V, f_i, post_ind, pre_ind, vel_ind, sign)

Kernel function for the computation of the Neumann-Computation for a
Non-Equilibrium Ansatz.
"""
# Define a functor for the bondary
function compute_neumann!{T<:Direction}(V::Neumann{T, NonEqBounce, D2Q9},
                         lattice::Lattice,
                        pre_ind::Array{Int64, 1},
                        post_ind::Array{Int64, 1},
                        vel_ind::Array{Int64, 1},
                         sign::Bool)

  f_i = lattice.f_temp
  # Compute the new values
  if sign
    rho_i = (sum(f_i[vel_ind]) + 2. *
             sum(f_i[pre_ind])) / (1. - V.vel)

    ru = rho_i * V.vel
    f_i[post_ind[1]] = f_i[pre_ind[1]] - (2./3.) * ru 
    f_i[post_ind[2]] = f_i[pre_ind[2]] - (1./6.) * ru +
      0.5 * (f_i[vel_ind[2]] - f_i[vel_ind[3]])
    f_i[post_ind[3]] = f_i[pre_ind[3]] - (1./6.) * ru +
      0.5 * (f_i[vel_ind[3]] - f_i[vel_ind[2]])
  
  else
    rho_i = (sum(f_i[vel_ind]) + 2. *
             sum(f_i[pre_ind])) / (1. + V.vel)

    ru = rho_i * V.vel
    f_i[post_ind[1]] = f_i[pre_ind[1]] +
      (2./3.) * ru 
    f_i[post_ind[2]] = f_i[pre_ind[2]] + (1./6.) * ru -
      0.5 * (f_i[vel_ind[2]] - f_i[vel_ind[3]])
    f_i[post_ind[3]] = f_i[pre_ind[3]] + (1./6.) * ru -
      0.5 * (f_i[vel_ind[3]] - f_i[vel_ind[2]])

    end

end

"""
    boundary(grid, bound::Neumann{T<:Direction, NonEqBounce, D2Q9})

Implementation of the interface function for the
computation of the neumann boundary condition. Dispath
over the direction.
"""
function boundary(grid::Grid_2D,
                  bound::Neumann{North, NonEqBounce,
                                 D2Q9}, d2q9::D2Q9)
  
    for row in bound.rows, col in bound.cols
      # Call the generic function
      compute_neumann(bound, grid.lattice[row, col],
                            velset.dict[East], velset.dict[West],
                            [1, 3, 5], false)
    end 
  
end

function boundary(grid::Grid_2D, bound::Neumann{South, NonEqBounce,
                                 D2Q9})
  
  for row in bound.rows, col in bound.cols
    # Call the generic function
      compute_neumann(bound, grid.lattice[row, col],
                            velset.dict[East], velset.dict[West],
                            [1, 3, 5], false)
  end

end

function boundary!(grid::Grid_2D,
                  bound::Neumann{West, NonEqBounce,
                                 D2Q9}, d2q9::D2Q9)
  
    for row in bound.rows, col in bound.cols
      # Call the generic function
      compute_neumann(bound, grid.lattice[row, col],
                            velset.dict[East], velset.dict[West],
                            [1, 3, 5], false)
    end

end

function boundary!(grid::Grid_2D,
                   bound::Neumann{East, NonEqBounce, D2Q9},
                   velset::_2D)
  
    for row in bound.rows, col in bound.cols
      # Call the generic function
      compute_neumann(bound, grid.lattice[row, col],
                            velset.dict[West], velset.dict[East],
                            [1, 3, 5], false)
    end

end
