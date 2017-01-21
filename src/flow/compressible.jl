#!/usr/bin/env julia

# === Equilibrium Functions
function compute_f_eq(lbm::LBM{_2D, Compressible,
                               Streaming, Collision})

    for k in 1:lbm.grid.directions
      lbm.grid.f_eq[:,:,k] = f_eq(F, w_q[k], lbm.grid.density,
                                  lbm.grid.velocity[:,:,1].^2 +
                                  lbm.grid.velocity[:,:,2].^2,
                                  c_dot_uv(lbm.grid.velocity,
                                           c_x[k], c_y[k]))
    end
end

function f_eq(w::Float64, rho::Float64,
                u_2::Array{Float64,1}, c_uv::Array{Float64,1})

     w .* (rho .+ 3.0 .* c_uv .+ 4.5 .* c_uv.^2 .- 1.5 .* u_2)
  end
  
function f_eq(w::Float64,
              rho::Array{Float64, 2},
              u_2::Array{Float64,2}, c_uv::Array{Float64,2})
  w .* (rho .+ 3.0 .* c_uv  + 4.5 .* c_uv.^2 - 1.5 .* u_2)
  end



