#!/usr/bin/env julia

# === Equilibrium Functions
function compute_f_eq(grid::Grid_2D, d2q9::D2Q9{Compressible})

    for k = 1:9
        grid.f_eq[:,:,k] = f_eq(d2q9.w[k], grid.density,
                                grid.velocity[:,:,1].^2 +
                                grid.velocity[:,:,2].^2,
                                c_dot_uv(grid.velocity,
                                          d2q9.c_x[k], d2q9.c_y[k]))
    end
end

function f_eq(w::Float64, rho::Float64,
              u_2::Array{Float64,1}, c_uv::Array{Float64,1})

    w .* (rho .+ 3.0 .* c_uv .+ 4.5 .* c_uv.^2 .- 1.5 .* u_2)
end

function f_eq(w::Float64,
              rho::Array{Float64, 2},
              u_2::Array{Float64,2}, c_uv::Array{Float64,2})
    w .* (rho .+ 3.0 .* c_uv  + 4.5 .* c_uv.^2 - 1.5 .* u_2)
end

"""
   compute(lbm, name, time_step [,write_inc=0])

Computes all steps defined by the Array `time_step' on the lattice boltzmann
model specified by `lbm'. If `write_inc' is 0, only the first and last
time_step are stored. Else every `write_inc' step is written into a 'vtr' file.

"""
function compute(grid::Grid_2D, d2q9::D2Q9{Compressible},
                 collision::Collision,
                 stream::Array{Streaming, 1}, 
                 bound::Array{Boundary, 1}, name::String,
                 time_step::Float64, write_inc::Int64=0)

  # Setting the directory
  src_dir = pwd()
  vtk_dir = string(src_dir, "/vtk")

  mkpath(vtk_dir)
  progress = Progress(length(time_step), 1, "Computing...", 30)
  
  println("INITIALISE LATTICE NODES...")
  init_lattice_state(grid, d2q9)

  i = 1
  write_vtk(grid, name, i)
  next!(progress)

    if write_inc == 0
        for i_step = 1.:time_step
            step(grid, d2q9, collision, stream, bound)

            next!(progress)
        end

        write_vtk(grid, name, 2)
    else
        cd(vtk_dir)
        for i_stel = 1.:time_step
            step(grid, d2q9, collision, stream, bound)

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
function step(grid::Grid_2D, d2q9::D2Q9, collision::Collision,
              stream::Array{Streaming, 1}, bound::Array{Boundary, 1})

    compute_collision(grid, collision)
    compute_streaming(grid, stream, d2q9)
    compute_boundary(grid, bound)
    compute_macro_var(grid, d2q9)
    compute_f_eq(grid, d2q9)

end

export compute

