#! /usr/bin/env julia
"""
This file contains the description of function relating to the vtr output file writing
"""

# =========== Output file writing

"""
    write_vtk(grid, name, step)

Writes the current status of the macroskopic variables
into a .vtr file. Julia uses column based storing of
the data. X & Y direction have to be swapped.
"""
function write_vtk(grid::Grid_2D, name::String, step::Int64)
    # Idea write the values of the particles in to arrays
    # and write them into the file.

    file_name = string(name, replace(string(step), ".","_"))

    vtk_f = vtk_grid(file_name, grid.x_point,
                     grid.y_point, append=false)

    vtk_point_data(vtk_f, grid.density, "Density")
    vtk_velocity = grid.velocity ./ grid.density

    # Velocities are swapped! Since julia uses column major formats..
    vtk_point_data(vtk_f, vtk_velocity[:, :, 1], "y-Velocity")
    vtk_point_data(vtk_f, vtk_velocity[:, :, 2], "x-Velocity")

    vtk_save(vtk_f)

end


# =========== Util functions ===========

function replace_brackets(str::String)
  
  return replace(replace(str, "]", " "), "[", " ")

end

function replace_comma(str::String)
  return replace(str, ",", " ")
end
