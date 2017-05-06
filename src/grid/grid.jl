#! /usr/bin/env julia

using .Abstract_LBM


function get_next_index(grid::Grid, i::Int64, dir::Int64, c_s_i::Int64)
    i_n = 1 + mod(i - 1 + c_s_i + size(grid.f_prop)[dir], size(grid.f_prop)[dir])

    return i_n
end

function get_axis_vec(width, length, consts::LBM_Constants)
  x = linspace(0, width, width) * consts.phys_x / width
  y = linspace(0, length, length) * consts.phys_y / length

  return x, y

end

include("population.jl")
include("macro_var_block.jl")
include("single_2d_grid.jl")

export
    SingleGrid2D
