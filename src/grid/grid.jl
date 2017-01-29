#! /usr/bin/env julia

using .Abstract_LBM

type Grid_2D <: Grid

    x_point::Array{Float64, 1}
    y_point::Array{Float64, 1}
    width::Int64
    length::Int64
    directions::Int64
    f_prop::Array{Float64, 3}
    f_eq::Array{Float64, 3}
    f_temp::Array{Float64, 3}
    density::Array{Float64,2}
    velocity::Array{Float64,3}

    Grid_2D(consts::LBM_Constants, width::Int64,
            length::Int64, directions::Int64) =
        (
            (x, y) = get_axis_vec(width, length, consts);
            f_prop = zeros(width, length, directions);
            density = zeros(width, length); # Init with ones
            velocity = zeros(width, length, 2);
            f_eq = zeros(width, length, directions);
            f_temp = zeros(width, length, directions);
            new(x, y, width, length, directions, f_prop,  
                f_eq, f_temp, density, velocity);
        )
end

function get_next_index(grid::Grid, i::Int64, dir::Int64, c_s_i::Int64)
    i_n = 1 + mod(i - 1 + c_s_i + size(grid.f_prop)[dir], size(grid.f_prop)[dir])

    return i_n
end

function get_axis_vec(width, length, consts::LBM_Constants)
  x = linspace(0, width, width) * consts.phys_x / width
  y = linspace(0, length, length) * consts.phys_y / length

  return x, y

end

"""
Compute the initial values of the grid. Gets called before the first normal
iteration.
"""
function init_lattice_state!(grid::Grid_2D, velset::_2D)

  # The Initial values for the grid 
  for k in 1:length(velset.w), i in 1:grid.width, j in 1:grid.length 
    grid.f_eq[i, j, k] = velset.w[k]
  end

  grid.f_temp = grid.f_eq
  grid.f_prop = grid.f_eq
  compute_macro_var!(grid, velset)

end

export Grid_2D
