#! /usr/bin/env julia

# Type definition for the lattice boltzmann 2D mesh
type LBM_Incompressible{Velocity_Set, Incompressible,
                        Streaming, Collision}<: LBM
    
    constants::LBM_Constants
    grid::Grid{Velocity_Set}
    bound::Array{Boundary, 1} # For Bounces etc...


end

function init_lattice_state(lbm::LBM_Incompressible{D2Q9, Incompressible,
                                                    FullPeriodicStreaming, Collision})


  #The Initial values for the grid
  for k in 1:lbm.grid.directions, i in 1:lbm.grid.width, j in 1:lbm.grid.length 
    lbm.grid.f_eq[i, j, k] = copy(_D2Q9.w[k])
  end

  lbm.grid.f_temp = copy(lbm.grid.f_eq)
  lbm.grid.f_prop = copy(lbm.grid.f_eq)
  compute_macro_var(lbm)
end

function compute_f_eq(lbm::LBM{_2D, Incompressible,
                               Streaming, Collision})

    for k in 1:lbm.grid.directions
        lbm.grid.f_eq[:,:,k] = f_eq(F, w_q[k], lbm.grid.density,
                                    c_dot_uv(lbm.grid.velocity,
                                             c_x[k], c_y[k]))
    end
end


function f_eq(w::Float64, rho::Float64,
                c_uv::Array{Float64,1})
      w .* (rho .+ 3.0 .* c_uv)
  end
  
function f_eq(w::Float64,
              rho::Array{Float64, 2}, c_uv::Array{Float64,2})
      w .* (rho .+ 3.0 .* c_uv)
  end

function compute_macro_var(lbm::LBM_Incompressible{D2Q9, Incompressible,
                                                   FullPeriodicStreaming, Collision})

    for i in 1:lbm.grid.width, j in 1:lbm.grid.length
        lbm.grid.density[i,j] = density(lbm.grid.f_prop[i,j, :])
        lbm.grid.velocity[i, j, :] =
            Array([sum(lbm.grid.f_prop[i,j, [2 6 9]]) - sum(lbm.grid.f_prop[i,j, [4 7 8]]),
                  sum(lbm.grid.f_prop[i,j, [3 6 7]]) - sum(lbm.grid.f_prop[i,j, [5 8 9]])])

    end

end
