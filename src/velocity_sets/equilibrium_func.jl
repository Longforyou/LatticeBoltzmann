#! /usr/bin/env julia

# """
# This file contains a description of usefull function
#  for the computation of the equilibrium function
# """
# =========== Equilibrium Distribution function

 @acc function c_dot_uv(velo::Array{Float64,2}, 
                  c_x::Float64, c_y::Float64)
  c_x .* velo[:,1] .+ c_y .* velo[:,2]
  end
  
  @acc function c_dot_uv(velo::Array{Float64,3}, 
                    c_x::Float64, c_y::Float64)
    c_x .* velo[:,:,1] .+ c_y .* velo[:,:,2]
  end

function compute_f_eq(grid::Grid_2D, d2q9::D2Q9{Incompressible})

    @parallel for k = 1:length(w)
        @inbounds @fastmath grid.f_eq[:,:,k] =
            f_eq(d2q9.w[k], grid.density,
                 c_dot_uv(grid.velocity,
                          d2q9.c_x[k], d2q9.c_y[k]))
    end
end

function f_eq(w::Float64, rho::Float64, c_uv::Array{Float64,1})
    w .* (rho .+ 3.0 .* c_uv)
end

function f_eq(w::Float64, rho::Array{Float64, 2}, c_uv::Array{Float64,2})
    w .* (rho .+ 3.0 .* c_uv)
end

# === Equilibrium Functions
function compute_f_eq(grid::Grid_2D, d2q9::D2Q9{Compressible})

    @parallel for k = 1:9
        @inbounds @fastmath grid.f_eq[:,:,k] = f_eq(d2q9.w[k], grid.density,
                                grid.velocity[:,:,1].^2 +
                                grid.velocity[:,:,2].^2,
                                c_dot_uv(grid.velocity,
                                          d2q9.c_x[k], d2q9.c_y[k]))
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
