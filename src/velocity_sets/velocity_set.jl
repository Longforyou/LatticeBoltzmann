#! /usr/bin/env julia

# using .Abstract_LBM

"""
This module contains definitions of types for referencing
the different kinds of velocity fields.
Example
-------
using velocity_set._D2Q9

make the weights w and the direction arrays c_x & c_y
 visible in the global namespace
"""

module _D2Q9
  using ..Abstract_LBM

  abstract D2Q9 <: _2D

  # Some constant fields
  const c_x = Array{Float64,1}([0., 1., 0., -1., 0., 1., -1., -1., 1.])
  const c_y = Array{Float64,1}([0., 0., 1., 0., -1., 1., 1., -1., -1.])
  const w = Array{Float64,1}([4/9, 1/9, 1/9, 1/9, 1/9,
                              1/36, 1/36, 1/36, 1/36])

  # ===== Macro Vars =================================
  function velo(V::D2Q9, f_prop::Array{Float64, 1})

      return Array([sum(f_prop[[2 6 9]]) - sum(f_prop[[4 7 8]]),
                  sum(f_prop[[3 6 7]]) - sum(f_prop[[5 8 9]])])
  end


  # Make the variables visible in the global namespace
  export c_x, c_y, w, D2Q9, velo

end # module _D2Q9

# Include the 2D velocity sets
# include("velocity_2d.jl"); export _D2Q9

# TODO 3d velocity sets

include("equilibrium_func.jl")
include("macro_var.jl")

