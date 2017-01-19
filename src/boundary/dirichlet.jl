#! /usr/bin/env julia


# ======= Pressure Dirichlet Conditions
abstract Dirichlet <: Boundary

type Pressure{T <: Direction, S <: SolType, V<:Velocity_Set} <: Boundary
  rho::Float64
  rows::Array{Int64, 1}
  cols::Array{Int64, 1}
   

end

function compute_dirichlet(V::D2Q9, f_i::Array{Float64, 1}, rho::Float64,
                           pre_ind::Array{Int64, 1}, 
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
function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::Pressure{North, NonEqBounce})
  
  for row in bound.rows, col in bound.cols
    # Call the generitc function
    lbm.grid.f_temp[row, col, :] =
        compute_dirichlet(lbm.grid.f_temp[row, col, :],
                          bound.rho, [4, 7, 8], [2, 5, 6],
                          [9, 1, 3], true )
  end 
  
end

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9}, 
                  bound::Pressure{South, NonEqBounce})
  
  for row in bound.rows, col in bound.cols
    # Call the generic function
      lbm.grid.f_temp[row, col, :] =
          compute_dirichlet(lbm.grid.f_temp[row, col, :],
                            bound.rho, [2, 5, 6], [4, 7, 8],
                            [9, 1, 3], true)
  end

end

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9}, 
                  bound::Pressure{West, NonEqBounce})


    for row in bound.rows, col in bound.cols
    # Call the generic function
      lbm.grid.f_temp[row, col, :] =
          compute_dirichlet(lbm.grid.f_temp[row, col, :],
                            bound.rho, [3, 7, 6], [1, 5, 8],
                            [9, 2, 4], false)
  end

end

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9}, 
                  bound::Pressure{East, NonEqBounce})

  for row in bound.rows, col in bound.cols
    # Call the generic function
    lbm.grid.f_temp[row, col, :] =
        compute_dirichlet(lbm.grid.f_temp[row, col, :],
                          bound.rho, [1, 5, 8], [3, 7, 6],
                          [9, 2, 4], true)
   
  end 
  
end
