#! /usr/bin/env julia

# This file contains a description of the grid, which contains a sparse octtree mesh data structure for storing the data.

mutable struct PatchGrid <: Grid
    patches::Array{Patch, 1}

end


mutable struct Grid_2D <: Grid

    x_point::Array{Float64, 1}
    y_point::Array{Float64, 1}
    width::Int64
    length::Int64
    directions::Int64
    lattices::Array{Lattice, 2}

    Grid_2D(consts::LBM_Constants, width::Int64,
            length::Int64, directions::Int64) =
        (
            (x, y) = get_axis_vec(width, length, consts);
            lattices = Array{Lattice}(width, length);
            
            for i = 1:width
              for j = 1:length
                lattices[i, j] = Lattice(directions, 2);
              end
            end;
            new(x, y, width, length, directions, lattices);
        )
end

function get_next_index(grid::Grid_2D, i::Int64, dir::Int64)
  
    i < dir ? i + 1 : 1
end

function get_prev_index(grid::Grid_2D, i::Int64, dir::Int64)
    i > 1 ? i - 1 : dir
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

  println("INITIALISE LATTICE NODES...")
  
  # The Initial values for the grid 
  for i in 1:grid.width, j in 1:grid.length 
      grid.lattices[i, j].f_eq = copy(velset.w)
      grid.lattices[i, j].f_temp = copy(velset.w)
      grid.lattices[i, j].f_prop = copy(velset.w)
  end

  compute_macro_var!(grid, velset)

end

include("macro_var_block.jl")
include("single_2d_grid.jl")

export
    SingleGrid2D
