#! /usr/bin/env julia


# ======= Pressure Dirichlet Conditions
abstract Dirichlet <: Boundary

type Pressure{T <: Direction, S <: SolType,
              V <: velocity_set} <: Boundary
  rho::Float64
  rows::Array{Int64, 1}
  cols::Array{Int64, 1}

end

function compute_dirichlet(V::velocity._D2Q9.D2Q9,
                           f_i::Array{Float64, 1},
                           rho::Float64, pre_ind::Array{Int64, 1}, 
                           post_ind::Array{Int64, 1},
                           rho_ind::Array{Int64, 1}, sign::Bool)

    ru = rho - abs(sum(f_i[rho_ind]) + 2. *
                   sum(f_i[pre_ind]))

    if sign
      f_i[post_ind[1]] = f_i[pre_ind[1]] - 
        (2./3.) * ru 
      f_i[post_ind[2]] = f_i[pre_ind[2]] - 
        (1./6.) * ru + 0.5 * (f_i[rho_ind[2]] - f_i[rho_ind[3]])
      f_i[post_ind[3]] = f_i[pre_ind[3]] - 
        (1./6.) * ru + 0.5 * (f_i[rho_ind[3]] - f_i[rho_ind[2]])
    else
      f_i[post_ind[1]] = f_i[pre_ind[1]] + 
        (2./3.) * ru 
      f_i[post_ind[2]] = f_i[pre_ind[2]] + 
        (1./6.) * ru + 0.5 * (f_i[rho_ind[2]] - f_i[rho_ind[3]])
      f_i[post_ind[3]] = f_i[pre_ind[3]] + 
        (1./6.) * ru + 0.5 * (f_i[rho_ind[3]] - f_i[rho_ind[2]])
    end

    return f_i
end

# ===========================================================
# === Dirichlet(Pressure) solution schemes
# ===========================================================
function boundary(lbm::LBM{V <: velocity_set._D2Q9.D2Q9, F <: Flow,
                           S <: Streaming, C <: Collision},
                  bound::Pressure{North, NonEqBounce,
                                  velocity_set._D2Q9.D2Q9})
  
  for row in bound.rows, col in bound.cols
    # Call the generitc function
    lbm.grid.f_temp[row, col, :] =
        compute_dirichlet(V, lbm.grid.f_temp[row, col, :],
                          bound.rho, [5, 8, 9], [3, 6, 7],
                          [1, 2, 4], true )
  end 
  
end

function boundary(lbm::LBM{V <: velocity_set._D2Q9.D2Q9, F <: Flow,
                           S <: Streaming, C <: Collision}, 
                  bound::Pressure{South, NonEqBounce,
                                  velocity_set._D2Q9.D2Q9})
  
  for row in bound.rows, col in bound.cols
      # Call the generic function
      lbm.grid.f_temp[row, col, :] =
          compute_dirichlet(V, lbm.grid.f_temp[row, col, :],
                            bound.rho, [3, 6, 7], [5, 8, 9],
                            [1, 2, 4], true)
  end

end

function boundary(lbm::LBM{V <: velocity_set._D2Q9.D2Q9, F <: Flow,
                           S <: Streaming, C <: Collision},
                  bound::Pressure{West, NonEqBounce,
                                  velocity_set._D2Q9.D2Q9})


    for row in bound.rows, col in bound.cols
    # Call the generic function
      lbm.grid.f_temp[row, col, :] =
          compute_dirichlet(V, lbm.grid.f_temp[row, col, :],
                            bound.rho, [4, 8, 7], [2, 6, 9],
                            [1, 3, 5], false)
  end

end

function boundary(lbm::LBM{V <: velocity_set._D2Q9.D2Q9, F <: Flow,
                           S <: Streaming, C <: Collision},
                  bound::Pressure{East, NonEqBounce,
                                  velocity_set._D2Q9.D2Q9})

  for row in bound.rows, col in bound.cols
      # Call the generic function
      lbm.grid.f_temp[row, col, :] =
          compute_dirichlet(V, lbm.grid.f_temp[row, col, :],
                              bound.rho, [2, 6, 9], [4, 8, 7],
                              [1, 3, 5], true)
   
  end 
  
end
