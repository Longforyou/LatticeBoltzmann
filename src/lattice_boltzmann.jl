#! /usr/bin/env julia

function compute_collision(lbm::Lattice_Boltzmann_2D{Cells.D2Q9})

  # The collision of the particles is indepentend of the particle type
  lbm.grid.f_temp = collision(lbm.grid.f_prop, lbm.constants.omega, 
                              lbm.grid.f_eq)

end

# =========== Propagations
function compute_propagation(lbm::Lattice_Boltzmann_2D{Cells.D2Q9})
    for bound in lbm.propagation
        boundary(lbm, bound)
    end
end


# =========== Steaming

function compute_streaming(lbm::Lattice_Boltzmann_2D{Cells.D2Q9})
    #c_x::Array{Float64, 1}, c_y::Array{Float64, 1})


  # Distribution direction
  lbm.grid.f_prop[:, :, 1] = circshift(lbm.grid.f_temp[:, : ,1],[ 0  1])
  lbm.grid.f_prop[:, :, 2] = circshift(lbm.grid.f_temp[:, :, 2],[ 1  0])
  lbm.grid.f_prop[:, :, 3] = circshift(lbm.grid.f_temp[:, :, 3],[ 0 -1])
  lbm.grid.f_prop[:, :, 4] = circshift(lbm.grid.f_temp[:, :, 4],[-1  0])
  lbm.grid.f_prop[:, :, 5] = circshift(lbm.grid.f_temp[:, :, 5],[ 1  1])
  lbm.grid.f_prop[:, :, 6] = circshift(lbm.grid.f_temp[:, :, 6],[ 1 -1])
  lbm.grid.f_prop[:, :, 7] = circshift(lbm.grid.f_temp[:, :, 7],[-1 -1])
  lbm.grid.f_prop[:, :, 8] = circshift(lbm.grid.f_temp[:, :, 8],[-1  1])
  lbm.grid.f_prop[:, :, 9] = lbm.grid.f_temp[:, :, 9]
  
end

# ===========
# ==== Boundary Functions
# ===========

function compute_boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9})

  for bound in lbm.bound
    boundary(lbm, bound)
  end
end

# Neumann solution schemes
function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::Neumann{North, NonEqBounce})
  
    for row in bound.rows, col in bound.cols
      # Call the generic function
      lbm.grid.f_temp[row, col, :] =
          compute_neumann(lbm.grid.f_temp[row, col, :],
                          bound.vel, [4, 7, 8], [2, 5, 6],
                          [9, 1, 3], true )
    end 

    # Correct the velocity on the boundary
    # lbm.grid.velocity[bound.rows, bound.cols, :] = bound.vel
  
end

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9}, 
                  bound::Neumann{South, NonEqBounce})
  
  for row in bound.rows, col in bound.cols
    # Call the generic function
    lbm.grid.f_temp[row, col, :] =
        compute_neumann(lbm.grid.f_temp[row, col, :],
                        bound.vel, [2, 5, 6], [4, 7, 8],
                        [9, 1, 3], false)
  end
 
  #Correct the velocity on the boundary
  # lbm.grid.velocity[bound.rows, bound.cols, :] = bound.vel
  

end

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9}, 
                  bound::Neumann{West, NonEqBounce})
  
    for row in bound.rows, col in bound.cols
      # Call the generic function
        lbm.grid.f_temp[row, col, :] =
            compute_neumann(lbm.grid.f_temp[row, col, :],
                            bound.vel, [3, 7, 6], [1, 5, 8],
                            [9, 2, 4], false)
    end
 
    # Correct the velocity on the boundary
    lbm.grid.velocity[bound.rows, bound.cols, :] = bound.vel
  

end

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9}, 
                  bound::Neumann{East, NonEqBounce})
 
  for row in bound.rows, col in bound.cols
    # Call the generic function
    lbm.grid.f_temp[row, col, :] =
        compute_neumann(lbm.grid.f_temp[row, col, :],
                        bound.vel, [1, 5, 8], [3, 7, 6],
                        [9, 2, 4], true)
   
  end 
  
  # Correct the velocity on the boundary
  # lbm.grid.velocity[bound.rows, bound.cols, :] = bound.vel

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
