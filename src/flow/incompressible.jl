#! /usr/bin/env julia

abstract Incompressible <: Flow

# Type definition for the lattice boltzmann 2D mesh
type LBM_Incompressible{V, F, S, C} <: LBM{V::Velocity_Set,
                                           F::Incompressible,
                                           S::Streaming,
                                           C::Collision}
  
  constants::LBM_Constants
  grid::Grid{V}
  bound::Array{Boundary, 1} # For Bounces etc...

  LBM_Incompressible(constants::LBM_Constants,
                       grid::Grid{V},
                       bound::Array{Boundary, 1}) =
    new(constants, grid, bound)

end

function compute_f_eq(lbm::LBM{V::velocity_set._D2Q9,
                               F::Incompressible,
                               S::Streaming, C::Collision})

    for k in 1:lbm.grid.directions
        lbm.grid.f_eq[:,:,k] = f_eq(F, w_q[k], lbm.grid.density,
                                    c_dot_uv(lbm.grid.velocity,
                                             c_x[k], c_y[k]))
    end
end

  function f_eq(F::Incompressible, w::Float64, rho::Float64, 
                c_uv::Array{Float64,1})
      w .* (rho .+ 3.0 .* c_uv)
  end
  
function f_eq(F::Incompressible, w::Float64,
              rho::Array{Float64, 2}, c_uv::Array{Float64,2})
      w .* (rho .+ 3.0 .* c_uv)
  end

