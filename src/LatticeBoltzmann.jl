#! /usr/bin/env julia

#__precompile__()

module LatticeBoltzmann

using  ProgressMeter, WriteVTK, ParallelAccelerator


# Process all other files
include("setup.jl")

"""
   compute(lbm, name, time_step [,write_inc=0])

Computes all steps defined by the Array `time_step' on the lattice boltzmann
model specified by `lbm'. If `write_inc' is 0, only the first and last
time_step are stored. Else every `write_inc' step is written into a 'vtr' file.

"""
function compute(grid::Grid, collision::Collision, stream::Streaming,
                 bound::Array{Boundary, 1}, name::String,
                 time_step::Array{Float64}, write_inc::Int64=0)

  # Setting the directory
  src_dir = pwd()
  vtk_dir = string(src_dir, "/vtk")

  mkpath(vtk_dir)
  progress = Progress(length(time_step), 1, "Computing...", 30)
  
  println("INITIALISE LATTICE NODES...")
  init_lattice_state(lbm)

  i = 1
  write_vtk(lbm, name, i)
  next!(progress)

    if write_inc == 0
        for i_step in time_step
            step(grid, collision, stream, bound)

            next!(progress)
        end

        write_vtk(grid, name, 2)
    else
        cd(vtk_dir)
        for i_stel in time_step
            step(lbm)

            if i % write_inc == 0
                write_vtk(grid, name, i)
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
function step(grid::Grid, collision::Collision,
              stream::Streaming, bound::Array{Boundary, 1})

    compute_collision(grid, collision)
    compute_streaming(grid, stream)
    compute_boundary(grid, bound)
    compute_macro_var(grid)
    compute_f_eq(grid)
>>>>>>> dev_newidea

end

export compute

end # module
