#! /usr/bin/env julia

# Type definition for the lattice boltzmann 2D mesh

function compute_f_eq(grid::Grid_2D, d2q9::D2Q9{Incompressible})

    for k = 1:length(w)
        grid.f_eq[:,:,k] = f_eq(d2q9.w[k], grid.density,
                                    c_dot_uv(grid.velocity,
                                             d2q9.c_x[k], d2q9.c_y[k]))
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
    init_lattice_state(lbm, w)

Compute the initial values of the grid. Gets called before the first normal
iteration.
"""
function init_lattice_state(grid::Grid_2D, d2q9::D2Q9)


  # The Initial values for the grid 
  for k = 1:9, i = 1:grid.width, j = 1:grid.length 
    grid.f_eq[i, j, k] = copy(d2q9.w[k])
  end

  grid.f_temp = copy(grid.f_eq)
  grid.f_prop = copy(grid.f_eq)
  compute_macro_var(grid, d2q9)


end

# ===========
