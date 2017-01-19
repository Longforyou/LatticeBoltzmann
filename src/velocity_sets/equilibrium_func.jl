#! /usr/bin/env julia

# =========== Equilibrium Distribution function

  function c_dot_uv(velo::Array{Float64,2}, 
                  c_x::Float64, c_y::Float64)
  c_x .* velo[:,1] .+ c_y .* velo[:,2]
  end
  
  function c_dot_uv(velo::Array{Float64,3}, 
                    c_x::Float64, c_y::Float64)
    c_x .* velo[:,:,1] .+ c_y .* velo[:,:,2]
  end

  
  function f_eq(w::Float64, rho::Float64, 
                u_2::Array{Float64,1}, c_uv::Array{Float64,1})
     w .* (rho .+ 3.0 .* c_uv .+ 4.5 .* c_uv.^2 .- 1.5 .* u_2)
  end
  
  function f_eq(w::Float64, rho::Array{Float64, 2},
                u_2::Array{Float64,2}, c_uv::Array{Float64,2})
  w .* (rho .+ 3.0 .* c_uv  + 4.5 .* c_uv.^2 - 1.5 .* u_2)
  end

  function compute_f_eq(lbm::LBM{D2Q9, F, S, C},
                         w_q::Array{Float64,1}, 
                         c_x::Array{Float64,1},
                         c_y::Array{Float64,1})

    for k in 1:lbm.grid.directions
      lbm.grid.f_eq[:,:,k] = f_eq(w_q[k], lbm.grid.density,
                                  lbm.grid.velocity[:,:,1].^2 +
                                  lbm.grid.velocity[:,:,2].^2,
                                  c_dot_uv(lbm.grid.velocity, c_x[k], c_y[k]))
    end
  end

# ===========
