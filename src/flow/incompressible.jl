#! /usr/bin/env julia

"""
Compute the initial values of the grid. Gets called before the first normal
iteration.
"""
function init_lattice_state(grid::Grid_2D)

  # The Initial values for the grid 
  for k in 1:length(w), i in 1:grid.width, j in 1:grid.length 
    grid.f_eq[i, j, k] = copy(w[k])
  end

  grid.f_temp = copy(grid.f_eq)
  grid.f_prop = copy(grid.f_eq)
  compute_macro_var(grid)


end
# ===========
