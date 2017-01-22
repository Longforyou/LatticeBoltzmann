#! /usr/bin/env julia

immutable PeriodicPressure{T <: Direction,
                      V <: _2D} <: Boundary

    rho_inlet::Float64
    rho_outlet::Float64
    inlet_row::Int64
    inlet_col::Array{Int64, 1}
    outlet_row::Int64
    outlet_col::Array{Int64, 1}

end

# ============================================================
# ==== Periodic Propagation
# ============================================================
function periodic_pressure(grid::Grid_2D, k::Int64,
                           bound_row::Int64,
                           bound_col::Array{Int64,1},
                           bound_rho::Float64, w::Float64,
                           c_x::Float64, c_y::Float64)

    f_eq(w, bound_rho,
         grid.velocity[bound_row, bound_col, 1].^2 .+
         grid.velocity[bound_row, bound_col, 2].^2, 
         c_dot_uv( grid.velocity[bound_row, bound_col, :], c_x, c_y)) .+ 
             ( grid.f_temp[bound_row, bound_col, k]
               .-  grid.f_eq[bound_row, bound_col, k])

end 

function boundary(grid::Grid_2D,
                  bound::PeriodicPressure{West, D2Q9})

    using _D2Q9: w, c_x, c_y

    # Compute the densities for the inlet and outlet,
    # where the pressure is known
    for k in 1:9
        #Inlet
        grid.f_temp[bound.inlet_row, bound.inlet_col, k] =
            periodic_pressure(grid, k, bound.outlet_row-1,
                              bound.outlet_col,
                              bound.rho_inlet, w[k], c_x[k], c_y[k])

        #Outlet
        grid.f_temp[bound.outlet_row, bound.outlet_col, k] =
            periodic_pressure(grid, k, bound.inlet_row+1,
                              bound.inlet_col,
                              bound.rho_outlet, w[k], c_x[k], c_y[k])

    end

end
