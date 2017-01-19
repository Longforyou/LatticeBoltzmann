#! /usr/bin/env julia
"""
This file contains the description of function relating to the vtr output file writing
"""
# =========== Output file writing

function write_vtk(lbm::Lattice_Boltzmann_2D, name::String, step::Int64)
  # Idea write the values of the particles in to arrays and write them into
  # the file.

  file_name = string(name, replace(string(step), ".","_"))
  x, y = get_axis_vec(lbm.grid, lbm.constants)

  vtk_f = vtk_grid(file_name, x, y, append=false)

  vtk_point_data(vtk_f, lbm.grid.density, "Density")
  vtk_velocity = lbm.grid.velocity ./ lbm.grid.density

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
