#! /usr/bin/env julia

using  ProgressMeter, WriteVTK

# Some constant fields
const c_x = Array{Float64,1}([1., 0., -1., 0., 1., -1., -1., 1., 0.])
const c_y = Array{Float64,1}([0., 1., 0., -1., 1., 1., -1., -1., 0.])
const w = Array{Float64,1}([1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36,
                                     1/36, 4/9])

# Type definition for the lattice boltzmann 2D mesh
type Lattice_Boltzmann_2D{T<:Cells.D2Q9}
  
  constants::LBM_Constants
  grid::Grid{T}
  bound::Array{Boundary, 1} # For Bounces etc...
  propagation::Array{Boundary, 1} # For Dirichlet Conditions

  Lattice_Boltzmann_2D(constants::LBM_Constants,
                       grid::Grid{T},
                       bound::Array{Boundary, 1}, prop::Array{Boundary, 1}) =
  new(constants, grid, bound, prop)

end

# ===========
"""
    init_lattice_state(lbm, w)

Compute the initial values of the grid. Gets called before the first normal
iteration.
"""
function init_lattice_state(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                            w::Array{Float64, 1})

  # The Initial values for the grid 
  for k in 1:lbm.grid.directions, i in 1:lbm.grid.width, j in 1:lbm.grid.length 
    lbm.grid.f_eq[i, j, k] = copy(w[k])
  end

  lbm.grid.f_temp = copy(lbm.grid.f_eq)
  lbm.grid.f_prop = copy(lbm.grid.f_eq)
  compute_macro_var(lbm)


end

# ===========
"""
   compute(lbm, name, time_step [,write_inc=0])

Computes all steps defined by the Array `time_step' on the lattice boltzmann
model specified by `lbm'. If `write_inc' is 0, only the first and last
time_step are stored. Else every `write_inc' step is written into a 'vtr' file.

"""
function compute(lbm::Lattice_Boltzmann_2D, name::String, 
                 time_step::Array{Float64}, write_inc::Int64=0)

  src_dir = pwd()
  vtk_dir = string(src_dir, "/vtk")

  mkpath(vtk_dir)
  progress = Progress(length(time_step), 1, "Computing...", 30)
  
  println("INITIALISE LATTICE NODES...")
  init_lattice_state(lbm, w)

  i = 1
  write_vtk(lbm, name, i)
  next!(progress)

    if write_inc == 0
        for i_step in time_step
            step(lbm)

            next!(progress)
        end

        write_vtk(lbm, name, 2)
    else
        cd(vtk_dir)
        for i_stel in time_step
            step(lbm)

            if i % write_inc == 0
                write_vtk(lbm, name, i)
            end

            i += 1
            next!(progress)
        end
        cd("..")
    end


end

# ===========
"""
    step(lbm)

Computes all individual steps needed in one iteration of the lattice boltzmann
scheme. 
  1. equilibrium distribution function
  2. collision operator ( only the BGK is implemented)
  3. periodic pressure conditions ( if there are any in the field `bound')
  4. steaming of the distributions to neighbouring nodes
  5. bounce backs on the values specified in all `BounceCondition' objects in
  `bound'
"""
function step(lbm::Lattice_Boltzmann_2D{Cells.D2Q9})

    compute_collision(lbm)
    compute_propagation(lbm)
    compute_streaming(lbm)

    compute_boundary(lbm)
    compute_macro_var(lbm)
    compute_f_eq(lbm, w, c_x, c_y)

end

# =========== Equilibrium Distribution function

  function c_dot_uv(velo::Array{Float64,2}, 
                  c_x::Float64, c_y::Float64)
  c_x .* velo[:,1] .+ c_y .* velo[:,2]
  end
  
  function c_dot_uv(velo::Array{Float64,3}, 
                    c_x::Float64, c_y::Float64)
    c_x .* velo[:,:,1] .+ c_y .* velo[:,:,2]
  end
  
  function f_eq(w::Float64, rho::Float64, 
                u_2::Array{Float64,1}, c_uv::Array{Float64,1})
     w .* (rho .+ 3.0 .* c_uv .+ 4.5 .* c_uv.^2 .- 1.5 .* u_2)
  end
  
  function f_eq(w::Float64, rho::Array{Float64, 2},
                u_2::Array{Float64,2}, c_uv::Array{Float64,2})
  w .* (rho .+ 3.0 .* c_uv  + 4.5 .* c_uv.^2 - 1.5 .* u_2)
  end

# write the function purely by passing arrays..
  function compute_f_eq(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                         w_q::Array{Float64,1}, 
                         c_x::Array{Float64,1},
                         c_y::Array{Float64,1})

    for k in 1:lbm.grid.directions
      lbm.grid.f_eq[:,:,k] = f_eq(w_q[k], lbm.grid.density,
                                  lbm.grid.velocity[:,:,1].^2 +
                                  lbm.grid.velocity[:,:,2].^2,
                                  c_dot_uv(lbm.grid.velocity, c_x[k], c_y[k]))
    end
  end

# ===========

function collision(f_prop::Array{Float64, 3}, omega::Float64,
                     f_eq::Array{Float64, 3}) 
    f_prop .* (1. - omega) .+ f_eq .* omega
end

function compute_collision(lbm::Lattice_Boltzmann_2D{Cells.D2Q9})

  # The collision of the particles is indepentend of the particle type
  lbm.grid.f_temp = collision(lbm.grid.f_prop, lbm.constants.omega, 
                              lbm.grid.f_eq)

end

# ===== Macro Vars ==========================================================
function velo(f_prop::Array{Float64, 1})

    return Array([sum(f_prop[[1 5 8]]) - sum(f_prop[[3 6 7]]),
                sum(f_prop[[2 5 6]]) - sum(f_prop[[4 7 8]])])
end

function density(f_prop::Array{Float64, 1})
    return sum(f_prop)
end

function compute_macro_var(lbm::Lattice_Boltzmann_2D{Cells.D2Q9})

      for i in 1:lbm.grid.width, j in 1:lbm.grid.length
          lbm.grid.density[i,j] = density(lbm.grid.f_prop[i,j, :])
          lbm.grid.velocity[i, j, :] = velo(lbm.grid.f_prop[i,j,:])
      end

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

# ============================================================
# ==== Corners
# ============================================================
function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::Corner)

  for row in bound.row, col in bound.col
    # Convenience pointer
    f_i = lbm.grid.f_temp[row, col, :]

    # Switch between the four cases
    if bound.quadrant == 0x01
        # f_3, f_7 & f_4 are known
        f_i[1] = f_i[3];  f_i[2] = f_i[4]
        f_i[5] = f_i[7];  f_i[6] = 0.; f_i[8] = 0.

    elseif bound.quadrant == 0x02
        # f_1, f_8 & f_4 are known
        f_i[3] = f_i[1];  f_i[2] = f_i[4]
        f_i[6] = f_i[8]; f_i[5] = 0.; f_i[7] = 0.

    elseif bound.quadrant == 0x03
        f_i[3] = f_i[1];  f_i[4] = f_i[2]
        f_i[7] = f_i[5];  f_i[5] = 0.; f_i[7] = 0.

    else
        f_i[1] = f_i[3];  f_i[4] = f_i[2]
        f_i[8] = f_i[6];  f_i[5] = 0.; f_i[7] = 0.

    end
    f_i[9] = bound.rho - sum(f_i[1:end-1])

    lbm.grid.f_temp[row, col, :] = f_i

   end
end

# ============================================================
# === Bounces
# ============================================================

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::Bounce{North})
    
    lbm.grid.f_prop[bound.rows, bound.cols, [4, 7, 8]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [2, 5, 6]]
end

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::Bounce{South})
    lbm.grid.f_prop[bound.rows, bound.cols, [2, 5, 6]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [4, 7, 8]]
end

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::Bounce{West})
    lbm.grid.f_prop[bound.rows, bound.cols, [3, 7, 6]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [1, 5, 8]]
end

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::Bounce{East})
    lbm.grid.f_prop[bound.rows, bound.cols, [1, 5, 8]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [3, 7, 6]]
end

# ============================================================
# === Open Boundary
# ============================================================

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::OpenBound{East})
    lbm.grid.f_prop[bound.rows, bound.cols, [1, 5, 8]] =
        2. * lbm.grid.f_prop[bound.rows-1, bound.cols, [1, 5, 8]] -
        lbm.grid.f_prop[bound.rows-2, bound.cols, [1, 5, 8]]

    lbm.grid.velocity[bound.rows, bound.cols, 2] = 0.
end

# ============================================================
# === Periodic Propagation
# ============================================================

function periodic_pressure(lbm::Lattice_Boltzmann_2D, k::Int64,
                                bound_row::Int64,
                                bound_col::Array{Int64,1},
                                bound_rho::Float64, w::Float64,
                                c_x::Float64, c_y::Float64)

  f_eq(w, bound_rho,
       lbm.grid.velocity[bound_row, bound_col, 1].^2 .+
           lbm.grid.velocity[bound_row, bound_col, 2].^2, 
  c_dot_uv(lbm.grid.velocity[bound_row, bound_col, :], c_x, c_y)) .+ 
  (lbm.grid.f_temp[bound_row, bound_col, k]
     .- lbm.grid.f_eq[bound_row, bound_col, k])

end 

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::PeriodicPressure{West})
  
  # Compute the densities for the inlet and outlet, where the pressure is
    # known
    for k in 1:lbm.grid.directions
      #Inlet
      lbm.grid.f_temp[bound.inlet_row, bound.inlet_col, k] =
        periodic_pressure(lbm, k, bound.outlet_row-1, bound.outlet_col,
                          bound.rho_inlet, w[k], c_x[k], c_y[k])

      #Outlet
      lbm.grid.f_temp[bound.outlet_row, bound.outlet_col, k] =
        periodic_pressure(lbm, k, bound.inlet_row+1, bound.inlet_col,
                          bound.rho_outlet, w[k], c_x[k], c_y[k])

    end

end

# =========== Output file writing

function write_vtk(lbm::Lattice_Boltzmann_2D, name::String, step::Int64)
  # Idea write the values of the particles in to arrays and write them into
  # the file.

  file_name = string(name, replace(string(step), ".","_"))
  x, y = get_axis_vec(lbm.grid, lbm.constants)

  vtk_f = vtk_grid(file_name, x, y, append=false)

  vtk_point_data(vtk_f, lbm.grid.density, "Density")
  vtk_velocity = lbm.grid.velocity ./ lbm.grid.density

  # Velocities are swapped! Since julia uses column major formats..
  vtk_point_data(vtk_f, vtk_velocity[:, :, 1], "y-Velocity")
  vtk_point_data(vtk_f, vtk_velocity[:, :, 2], "x-Velocity")


  vtk_save(vtk_f)

end


# =========== Util functions ====

function replace_brackets(str::String)
  
  return replace(replace(str, "]", " "), "[", " ")

end

function replace_comma(str::String)
  return replace(str, ",", " ")
end
