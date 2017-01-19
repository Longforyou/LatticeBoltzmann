#! /usr/bin/env julia

type PeriodicPressure{T <: Direction} <: Boundary
    rho_inlet::Float64
    rho_outlet::Float64
    inlet_row::Int64
    inlet_col::Array{Int64, 1}
    outlet_row::Int64
    outlet_col::Array{Int64, 1}

end

# ============================================================
# === Periodic Propagation
# ============================================================

function periodic_pressure(lbm::Lattice_Boltzmann_2D, k::Int64,
                                bound_row::Int64,
                                bound_col::Array{Int64,1},
                                bound_rho::Float64, w::Float64,
                                c_x::Float64, c_y::Float64)

  f_eq(w, bound_rho,
       lbm.grid.velocity[bound_row, bound_col, 1].^2 .+
           lbm.grid.velocity[bound_row, bound_col, 2].^2, 
  c_dot_uv(lbm.grid.velocity[bound_row, bound_col, :], c_x, c_y)) .+ 
  (lbm.grid.f_temp[bound_row, bound_col, k]
     .- lbm.grid.f_eq[bound_row, bound_col, k])

end 

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::PeriodicPressure{West})
  
  # Compute the densities for the inlet and outlet, where the pressure is
    # known
    for k in 1:lbm.grid.directions
      #Inlet
      lbm.grid.f_temp[bound.inlet_row, bound.inlet_col, k] =
        periodic_pressure(lbm, k, bound.outlet_row-1, bound.outlet_col,
                          bound.rho_inlet, w[k], c_x[k], c_y[k])

      #Outlet
      lbm.grid.f_temp[bound.outlet_row, bound.outlet_col, k] =
        periodic_pressure(lbm, k, bound.inlet_row+1, bound.inlet_col,
                          bound.rho_outlet, w[k], c_x[k], c_y[k])

    end

end
