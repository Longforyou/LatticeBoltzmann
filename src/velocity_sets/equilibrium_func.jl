#! /usr/bin/env julia

"""
This file contains a description of usefull function
 for the computation of the equilibrium function
"""
# =========== Equilibrium Distribution function

  function c_dot_uv(velo::Array{Float64,2}, 
                  c_x::Float64, c_y::Float64)
  c_x .* velo[:,1] .+ c_y .* velo[:,2]
  end
  
  function c_dot_uv(velo::Array{Float64,3}, 
                    c_x::Float64, c_y::Float64)
    c_x .* velo[:,:,1] .+ c_y .* velo[:,:,2]
  end

  
