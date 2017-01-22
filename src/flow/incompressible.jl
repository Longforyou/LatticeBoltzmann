#! /usr/bin/env julia

# Type definition for the lattice boltzmann 2D mesh

function compute_f_eq(grid::Grid_2D{D2Q9, Incompressible})

    for k in 1:length(w)
        grid.f_eq[:,:,k] = f_eq(_D2Q0.w[k], grid.density,
                                    c_dot_uv(grid.velocity,
                                             _D2Q0.c_x[k], _D2Q0.c_y[k]))
    end
end


function f_eq(w::Float64, rho::Float64, c_uv::Array{Float64,1})
    w .* (rho .+ 3.0 .* c_uv)
end

function f_eq(w::Float64, rho::Array{Float64, 2}, c_uv::Array{Float64,2})
    w .* (rho .+ 3.0 .* c_uv)
end


# ===========
"""
    init_lattice_state( w)

Compute the initial values of the grid. Gets called before the first normal
iteration.
"""
function init_lattice_state(grid::Grid_2D{D2Q9, Incompressible})

  using _D2Q9: w

  # The Initial values for the grid 
  for k in 1:length(w), i in 1:grid.width, j in 1:grid.length 
    grid.f_eq[i, j, k] = copy(w[k])
  end

  grid.f_temp = copy(grid.f_eq)
  grid.f_prop = copy(grid.f_eq)
  compute_macro_var(grid)


end

# ===========
