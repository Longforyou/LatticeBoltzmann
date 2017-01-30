#! /usr/bin/env julia

using ._D2Q9

# ====== Neumann conditions
# The generic type T, allows for different equilibirum implementations
immutable Neumann{T<:Direction, S<:SolType,
                  V <: Velocity_Set} <: Boundary
  vel::Float64
  rows::Array{Int64, 1}
  cols::Array{Int64, 1}

end

# Define a functor for the bondary
function compute_neumann{T<:Direction}(V::Neumann{T, NonEqBounce, D2Q9},
                         f_i::Array{Float64, 1},
                        post_ind::Array{Int64, 1},
                        pre_ind::Array{Int64, 1},
                        vel_ind::Array{Int64, 1},
                         sign::Bool)

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

  return f_i
end

# Neumann solution schemes for D2Q9
function boundary(grid::Grid_2D,
                  bound::Neumann{North, NonEqBounce,
                                 D2Q9}, d2q9::D2Q9)
  
    for row in bound.rows, col in bound.cols
      # Call the generic function
      grid.f_prop[row, col, :] =
            compute_neumann(bound, grid.f_temp[row, col, :],
                          bound.vel, velset.dict[South], velset.dict[North],
                          [1, 2, 4], true )
    end 
  
end

function boundary(grid::Grid_2D, bound::Neumann{South, NonEqBounce,
                                 D2Q9})
  
  for row in bound.rows, col in bound.cols
    # Call the generic function
    grid.f_prop[row, col, :] =
            compute_neumann(bound, grid.f_temp[row, col, :],
                        bound.vel, velset.dict[North], velset.dict[South],
                        [1, 2, 4], false)
  end

end

function boundary!(grid::Grid_2D,
                  bound::Neumann{West, NonEqBounce,
                                 D2Q9}, d2q9::D2Q9)
  
    for row in bound.rows, col in bound.cols
      # Call the generic function
        grid.f_prop[row, col, :] =
            compute_neumann(bound, grid.f_temp[row, col, :],
                            velset.dict[East], velset.dict[West],
                            [1, 3, 5], false)
    end
                            #bound.vel, [4, 8, 6], [2, 6, 9],

end

function boundary!(grid::Grid_2D,
                  bound::Neumann{East, NonEqBounce,
                                 D2Q9}, d2q9::D2Q9)
  
    for row in bound.rows, col in bound.cols
      # Call the generic function
        grid.f_prop[row, col, :] =
            compute_neumann(bound, grid.f_temp[row, col, :],
                            velset.dict[West], velset.dict[East],
                            [1, 3, 5], false)
    end

end
