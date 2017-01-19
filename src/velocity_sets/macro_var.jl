#! /usr/bin/env julia

# ===== Macro Vars ==========================================================
function velo(V::velocity_set._2D.D2Q9, f_prop::Array{Float64, 1})

    return Array([sum(f_prop[[2 6 9]]) - sum(f_prop[[4 7 8]]),
                sum(f_prop[[3 6 7]]) - sum(f_prop[[5 8 9]])])
end

function density(f_prop::Array{Float64, 1})
    return sum(f_prop)
end

function compute_macro_var(lbm::LBM{V <: Velocity_Set, F <: Flow,
                                    S <: Streaming, C <: Collision})

      for i in 1:lbm.grid.width, j in 1:lbm.grid.length
          lbm.grid.density[i,j] = density(lbm.grid.f_prop[i,j, :])
          lbm.grid.velocity[i, j, :] = velo(V, lbm.grid.f_prop[i,j,:])
      end

  end
