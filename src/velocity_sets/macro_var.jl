#! /usr/bin/env julia

# ===== Macro Vars ==========================================================
function velo(f_prop::Array{Float64, 1})

    return Array([sum(f_prop[[1 5 8]]) - sum(f_prop[[3 6 7]]),
                sum(f_prop[[2 5 6]]) - sum(f_prop[[4 7 8]])])
end

function density(f_prop::Array{Float64, 1})
    return sum(f_prop)
end

function compute_macro_var(lbm::Lattice_Boltzmann_2D{Cells.D2Q9})

      for i in 1:lbm.grid.width, j in 1:lbm.grid.length
          lbm.grid.density[i,j] = density(lbm.grid.f_prop[i,j, :])
          lbm.grid.velocity[i, j, :] = velo(lbm.grid.f_prop[i,j,:])
      end

  end
