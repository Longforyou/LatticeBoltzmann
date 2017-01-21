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
function (V::Neumann{Direction,
                     NonEqBounce, D2Q9})(f_i::Array{Float64, 1},
                        vel_ind::Array{Int64, 1},
                        pre_ind::Array{Int64, 1},
                        post_ind::Array{Int64, 1},
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
function boundary(lbm::LBM{D2Q9, Flow,
                           Streaming, Collision},
                  bound::Neumann{North, NonEqBounce,
                                 D2Q9})
  
    for row in bound.rows, col in bound.cols
      # Call the generic function
      lbm.grid.f_temp[row, col, :] =
          bound(lbm.grid.f_temp[row, col, :],
                          bound.vel, [5, 8, 9], [3, 6, 7],
                          [1, 2, 4], true )
    end 

    # Correct the velocity on the boundary
    # lbm.grid.velocity[bound.rows, bound.cols, :] = bound.vel
  
end

function boundary(lbm::LBM{D2Q9, Flow,
                           Streaming, Collision},
                  bound::Neumann{South, NonEqBounce,
                                 D2Q9})
  
  for row in bound.rows, col in bound.cols
    # Call the generic function
    lbm.grid.f_temp[row, col, :] =
        bound(lbm.grid.f_temp[row, col, :],
                        bound.vel, [3, 6, 7], [5, 8, 9],
                        [1, 2, 4], false)
  end
 
  #Correct the velocity on the boundary
  # lbm.grid.velocity[bound.rows, bound.cols, :] = bound.vel
  

end

function boundary(lbm::LBM{D2Q9, Flow,
                           Streaming, Collision},
                  bound::Neumann{West, NonEqBounce,
                                 D2Q9})
  
    for row in bound.rows, col in bound.cols
      # Call the generic function
        lbm.grid.f_temp[row, col, :] =
            bound(lbm.grid.f_temp[row, col, :],
                            bound.vel, [4, 8, 6], [2, 6, 9],
                            [1, 3, 5], false)
    end
 
    # Correct the velocity on the boundary
    lbm.grid.velocity[bound.rows, bound.cols, :] = bound.vel
  

end

function boundary(lbm::LBM{D2Q9, Flow,
                           Streaming, Collision},
                  bound::Neumann{Direction, NonEqBounce,
                                 D2Q9})
 
  for row in bound.rows, col in bound.cols
    # Call the generic function
    lbm.grid.f_temp[row, col, :] =
        bound(lbm.grid.f_temp[row, col, :], true)
   
  end 
  
  # Correct the velocity on the boundary
  # lbm.grid.velocity[bound.rows, bound.cols, :] = bound.vel

end

