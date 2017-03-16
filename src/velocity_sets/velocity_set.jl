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

immutable D2Q9{F<:Flow} <: _2D
    w::Array{Float64, 1}
    dict::Dict{DataType, Array{Int64, 1}}

    D2Q9() = (
      # Some constant fields
      new(Array{Float64, 1}([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]),
          Dict{DataType, Array{Int64, 1}}([(North, [3, 6, 7]),
           (South, [5, 8, 9]), (West, [4, 7, 8]), (East, [2, 6, 9])]))
    )

  end


function velo_1(f_prop::Array{Float64, 1})

    return f_prop[2] + f_prop[6] + f_prop[9] -
        f_prop[4] - f_prop[7] - f_prop[8]

end

function velo_2(f_prop::Array{Float64, 1})
    return f_prop[3] + f_prop[6] + f_prop[7] -
        f_prop[5] - f_prop[8] - f_prop[9]
end

  # Make the variables visible in the global namespace
  export D2Q9

end # module _D2Q9


using ._D2Q9

include("equilibrium_func.jl")
include("macro_var.jl")

