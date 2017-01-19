#! /usr/bin/env julia

__precompile__()

module LatticeBoltzmann

using  ProgressMeter, WriteVTK

# Process the abstract modelling
include("abstract_lbm.jl")

# Process all other files
include("setup.jl")

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
function compute(lbm::LBM, name::String, 
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
function step(lbm::LBM)

    compute_collision(lbm)
    compute_propagation(lbm)
    compute_streaming(lbm)

    compute_boundary(lbm)
    compute_macro_var(lbm)
    compute_f_eq(lbm, w, c_x, c_y)

end

export Lattice_Boltzmann_2D, compute

end # module
