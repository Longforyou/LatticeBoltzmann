#! /usr/bin/env julia

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

    vtk_density, vtk_x_vel, vtk_y_vel = get_macro_var(grid.lattices)
    
    vtk_point_data(vtk_f, vtk_density, "Density")
    vtk_point_data(vtk_f, vtk_x_vel, "x-Velocity")
    vtk_point_data(vtk_f, vtk_y_vel, "y-Velocity")

    vtk_save(vtk_f)

    # Free up the space
    vtk_density = 0
    vtk_x_vel = 0
    vtk_y_vel = 0

end

# =========== Util functions ===========

function replace_brackets(str::String)
  
  return replace(replace(str, "]", " "), "[", " ")

end

function replace_comma(str::String)
  return replace(str, ",", " ")
end
