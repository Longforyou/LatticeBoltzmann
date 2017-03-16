#! /usr/bin/env julia

__precompile__()

module LatticeBoltzmann

using  ProgressMeter, WriteVTK #, ParallelAccelerator


# Process all other files
include("setup.jl")

"""
   compute(grid, velset, collision, stream, bound, name, tim_step [,write_inc=0])

Computes all steps defined by the Array `time_step' on the lattice boltzmann
model specified by the following variables.
* grid::Grid | contains all Arrays for the computation
* velset::Velocity_Set | contains the description of the Velocity Setting
* collision::Collision | specifies the collision operator
* stream::Array{Streaming, 1} | contains the description of streaming operations
* bound::Array{Boundary, 1} | contains all boundary conditions of the models
 If `write_inc' is 0, only the first and last
time_step are stored. Else every `write_inc' step is written into a 'vtr' file. Not that this has an direct influence on the performance.

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
  
  println("INITIALISE LATTICE NODES...")
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
        @showprogress 3 "Computing..." for i_stel = 1.:time_step
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
   compute(grid, velset, collision, stream, bound, name, time_step
   analytical_sol, [,write_inc=0])

Computes all steps defined by the Array `time_step' on the lattice boltzmann
model specified by `lbm'. If `write_inc' is 0, only the first and last
time_step are stored. Else every `write_inc' step is written into a 'vtr' file.
This function does a error analysis for the passed analytical solution.

"""
function compute!(grid::Grid, velset::Velocity_Set,
                 collision::Collision,
                 stream::Array{Streaming, 1}, 
                 bound::Array{Boundary, 1}, name::String,
                 time_step::Float64, analytical_sol::Array{Float64, 1}, vel_dir::Int64,  write_inc::Int64=0)

  # Setting the directory
  src_dir = pwd()
  vtk_dir = string(src_dir, "/vtk")

    mkpath(vtk_dir)
    println("Run a LatticeBoltzmann-Simulation with ",
            time_step, " step.")

  # Generate the array for the convergence
  println(size(analytical_sol), size(grid.velocity[:, :, vel_dir]))
  conv_vec = zeros(Int64(time_step))
  
  println("INITIALISE LATTICE NODES...")
  init_lattice_state!(grid, velset)

  i = 1
  write_vtk(grid, name, i)

    if write_inc == 0
        @showprogress 3 "Computing..." for i_step = 1.:time_step
            step!(grid, velset, collision, stream, bound)
            conv_vec[Int64(i_step)] = norm((grid.velocity[:, :, vel_dir] ./ grid.density)' .- analytical_sol)

        end

        write_vtk(grid, name, 2)
    else
        cd(vtk_dir)
        @showprogress 3 "Computing..." for i_stel = 1.:time_step
            step!(grid, velset, collision, stream, bound)
            conv_vec[Int64(i_step)] = norm((grid.velocity[:, :, vel_dir] ./ grid.density)' .- analytical_sol)

            if i % write_inc == 0
                write_vtk(grid, name, i)
            end

            i += 1
        end
        cd("..")
    end

    return conv_vec

end

"""
    step(lbm)

Computes all individual steps needed in one iteration of the lattice boltzmann
scheme. 
  1. collision operator ( only the BGK is implemented)
  2. steaming of the distributions to neighbouring nodes
  3. computation of the passed boundary conditions
  4. computation of the macroskopic variables
  5. equilibrium distribution function
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
