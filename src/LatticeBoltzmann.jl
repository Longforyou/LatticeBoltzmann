#! /usr/bin/env julia

__precompile__()

module LatticeBoltzmann

using  ProgressMeter, WriteVTK #, ParallelAccelerator


# Process all other files
include("setup.jl")

"""
   compute(lbm, name, tim_step [,write_inc=0])

Computes all steps defined by the Array `time_step' on the lattice boltzmann
model specified by `lbm'. If `write_inc' is 0, only the first and last
time_step are stored. Else every `write_inc' step is written into a 'vtr' file.

"""
function compute!(grid::Grid, velset::Velocity_Set,
                 collision::Collision,
                 stream::Array{Streaming, 1}, 
                 bound::Array{Boundary, 1}, name::String,
                 time_step::Float64, write_inc::Int64=0)

  # Setting the directory
  src_dir = pwd()
  vtk_dir = string(src_dir, "/vtk")

    mkpath(vtk_dir)
    println("Run a LatticeBoltzmann-Simulation with ",
            time_step, " step.")
  
  init_lattice_state!(grid, velset)

  i = 1
  write_vtk(grid, name, i)

    if write_inc == 0
        @showprogress 3 "Computing..." for i_step = 1.:time_step
            step!(grid, velset, collision, stream, bound)

        end

        write_vtk(grid, name, 2)
    else
        cd(vtk_dir)
        @showprogress 3 "Computing..." for i_step = 1.:time_step
            step!(grid, velset, collision, stream, bound)

            if i % write_inc == 0
                write_vtk(grid, name, i)
            end

            i += 1
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
function step!(grid::Grid, velset::Velocity_Set, collision::Collision,
              stream::Array{Streaming, 1}, bound::Array{Boundary, 1})

    compute_collision!(grid, collision)
    compute_streaming!(grid, stream, velset)
    compute_boundary!(grid, bound, velset)
    compute_macro_var!(grid, velset)
    compute_f_eq!(grid, velset)

end

export compute!

end # module
