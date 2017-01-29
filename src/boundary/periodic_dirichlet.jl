#! /usr/bin/env julia
immutable PressurePeriodic_2D{T <: Direction,
                           V <: _2D} <: Boundary

    rho_inlet::Float64
    rho_outlet::Float64
    inlet_row::Int64
    inlet_col::Array{Int64, 1}
    outlet_row::Int64
    outlet_col::Array{Int64, 1}

end

# ============================================================
# ==== Periodic STreaming with Pressure condition
# ============================================================

function boundary!(grid::Grid_2D,
                   PFPS::PressurePeriodic_2D{West, D2Q9}, d2q9::D2Q9)

    # Compute the densities for the inlet and outlet,
    # where the pressure is known

    #Inlet
    @fastmath @inbounds grid.f_prop[PFPS.inlet_row, PFPS.inlet_col, :] =
        periodic_pressure(grid, d2q9, PFPS.outlet_row-1,
                          PFPS.outlet_col,
                          PFPS.rho_inlet)

    #Outlet
    @fastmath @inbounds grid.f_prop[PFPS.outlet_row, PFPS.outlet_col, :] =
        periodic_pressure(grid, d2q9, PFPS.inlet_row+1,
                          PFPS.inlet_col,
                          PFPS.rho_outlet)

end
