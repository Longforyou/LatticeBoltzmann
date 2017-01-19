#! /usr/bin/env julia


# ====== Neumann conditions
# The generic type T, allows for different equilibirum implementations
type Neumann{T<:Direction, S<:SolType, V <: Velocity_Set} <: Boundary
  vel::Float64
  rows::Array{Int64, 1}
  cols::Array{Int64, 1}

end

function compute_neumann(V::velocity_set._D2Q9.D2Q9,
                         f_i::Array{Float64, 1},
                         vel::Float64,
                         pre_ind::Array{Int64, 1},
                         post_ind::Array{Int64, 1},
                         vel_ind1::Array{Int64, 1}, sign::Bool)

  # Compute the new values
  if sign
    rho_i = (sum(f_i[vel_ind1]) + 2. *
             sum(f_i[pre_ind])) / (1. - vel)

    ru = rho_i * vel

    f_i[post_ind[1]] = f_i[pre_ind[1]] - (2./3.) * ru 
    f_i[post_ind[2]] = f_i[pre_ind[2]] - (1./6.) * ru +
      0.5 * (f_i[vel_ind1[2]] - f_i[vel_ind1[3]])
    f_i[post_ind[3]] = f_i[pre_ind[3]] - (1./6.) * ru +
      0.5 * (f_i[vel_ind1[3]] - f_i[vel_ind1[2]])
  
  else
    rho_i = (sum(f_i[vel_ind1]) + 2. *
             sum(f_i[pre_ind])) / (1. + vel)

    ru = rho_i * vel
    f_i[post_ind[1]] = f_i[pre_ind[1]] +
      (2./3.) * ru 
    f_i[post_ind[2]] = f_i[pre_ind[2]] + (1./6.) * ru -
      0.5 * (f_i[vel_ind1[2]] - f_i[vel_ind1[3]])
    f_i[post_ind[3]] = f_i[pre_ind[3]] + (1./6.) * ru -
      0.5 * (f_i[vel_ind1[3]] - f_i[vel_ind1[2]])

    end

  return f_i
end

# Neumann solution schemes for D2Q9
function boundary(lbm::LBM{V <: velocity._2D.D2Q9, F <: Flow,
                           S <: Streaming, C <: Collision},
                  bound::Neumann{North, NonEqBounce,
                                 velocity_set._2D.D2Q9})
  
    for row in bound.rows, col in bound.cols
      # Call the generic function
      lbm.grid.f_temp[row, col, :] =
          compute_neumann(V, lbm.grid.f_temp[row, col, :],
                          bound.vel, [5, 8, 9], [3, 6, 7],
                          [1, 2, 4], true )
    end 

    # Correct the velocity on the boundary
    # lbm.grid.velocity[bound.rows, bound.cols, :] = bound.vel
  
end

function boundary(lbm::LBM{V <: velocity._2D.D2Q9, F <: Flow,
                           S <: Streaming, C <: Collision},
                  bound::Neumann{South, NonEqBounce,
                                 velocity_set._2D.D2Q9})
  
  for row in bound.rows, col in bound.cols
    # Call the generic function
    lbm.grid.f_temp[row, col, :] =
        compute_neumann(V, lbm.grid.f_temp[row, col, :],
                        bound.vel, [3, 6, 7], [5, 8, 9],
                        [1, 2, 4], false)
  end
 
  #Correct the velocity on the boundary
  # lbm.grid.velocity[bound.rows, bound.cols, :] = bound.vel
  

end

function boundary(lbm::LBM{V <: velocity._2D.D2Q9, F <: Flow, S <: Streaming, C <: Collision}, 
                  bound::Neumann{West, NonEqBounce,
                                 velocity_set._2D.D2Q9})
  
    for row in bound.rows, col in bound.cols
      # Call the generic function
        lbm.grid.f_temp[row, col, :] =
            compute_neumann(V, lbm.grid.f_temp[row, col, :],
                            bound.vel, [4, 8, 6], [2, 6, 9],
                            [1, 3, 5], false)
    end
 
    # Correct the velocity on the boundary
    lbm.grid.velocity[bound.rows, bound.cols, :] = bound.vel
  

end

function boundary(lbm::LBM{V <: velocity._2D.D2Q9, F <: Flow,
                           S <: Streaming, C <: Collision},
                  bound::Neumann{East, NonEqBounce,
                                 velocity_set._2D.D2Q9})
 
  for row in bound.rows, col in bound.cols
    # Call the generic function
    lbm.grid.f_temp[row, col, :] =
        compute_neumann(V, lbm.grid.f_temp[row, col, :],
                        bound.vel, [2, 6, 9], [4, 8, 7],
                        [1, 3, 5], true)
   
  end 
  
  # Correct the velocity on the boundary
  # lbm.grid.velocity[bound.rows, bound.cols, :] = bound.vel

end

