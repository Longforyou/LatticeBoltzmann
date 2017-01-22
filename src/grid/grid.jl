#! /usr/bin/env julia

using .Abstract_LBM

type Grid_2D{V <: Velocity_Set}

    x_point::Array{Float64, 1}
    y_point::Array{Float64, 1}
    width::Int64
    length::Int64
    directions::Int64
    f_prop::Array{Float64, 3}
    f_eq::Array{Float64, 3}
    density::Array{Float64,2}
    velocity::Array{Float64,3}
    f_temp::Array{Float64, 3}

    Grid_2D(consts, width, length, directions) =
        (
            x, y = get_axis_vec(width, length, consts);
            f_prop = zeros(width, length, directions);
            density = zeros(width, length); # Init with ones
            velocity = zeros(width, length, 2);
            f_eq = zeros(width, length, directions);
            f_temp = zeros(width, length, directions);

            new(x, y, width, length, directions, f_prop,  
                f_eq, density, velocity, f_temp);
        )
end

function getNeighbours(grid::Grid, i::Int, j::Int)
    #=
    Returns the indices for referencing the 8 Lattice neighbours.
    =#

    neighbours = Array{Array{Int, 1},1}(8)

    i_p, i_n = get_dir_1_indices(grid, i)
    j_p, j_n = get_dir_2_indices(grid, j)

    neighbours[1] = [i_n, j]
    neighbours[2] = [i, j_p]
    neighbours[3] = [i_p, j]
    neighbours[4] = [i, j_n]
    neighbours[5] = [i_n, j_p]
    neighbours[6] = [i_p, j_p]
    neighbours[7] = [i_p, j_n]
    neighbours[8] = [i_n, j_n]

    return neighbours
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
