#! /usr/bin/env julia

immutable FullPeriodicStreaming_2D <: Streaming
    rows::UnitRange{Int64}
    cols::UnitRange{Int64}

    function FullPeriodicStreaming_2D(grid::Grid_2D)
        new(UnitRange{Int64}(1:grid.width),
            UnitRange{Int64}(1:grid.length))
    end

end

abstract InnerStreaming <: Streaming

# =========== Streaming
function compute_streaming!(grid::Grid, stream::Array{Streaming, 1}, velset::Velocity_Set)
    for stre in stream
        # println("Pre ", stre, "\nPopulation\n", grid.f_prop)
        @inbounds streaming!(stre, grid, velset)
        # println("Post ", stre, "\nPopulation\n", grid.f_prop)
    end
end

function streaming!(FPS::FullPeriodicStreaming_2D, grid::Grid_2D,
                    d2q9::D2Q9)

    # Distribution direction
    # grid.f_prop[FPS.rows, FPS.cols, 1] = grid.f_temp[FPS.rows, FPS.cols, 1]
    # grid.f_prop[FPS.rows, FPS.cols, 2] = circshift(grid.f_temp[FPS.rows, FPS.cols, 2],[ 0  1])
    # grid.f_prop[FPS.rows, FPS.cols, 3] = circshift(grid.f_temp[FPS.rows, FPS.cols, 3],[ 1  0])
    # grid.f_prop[FPS.rows, FPS.cols, 4] = circshift(grid.f_temp[FPS.rows, FPS.cols, 4],[ 0 -1])
    # grid.f_prop[FPS.rows, FPS.cols, 5] = circshift(grid.f_temp[FPS.rows, FPS.cols, 5],[-1  0])
    # grid.f_prop[FPS.rows, FPS.cols, 6] = circshift(grid.f_temp[FPS.rows, FPS.cols, 6],[ 1  1])
    # grid.f_prop[FPS.rows, FPS.cols, 7] = circshift(grid.f_temp[FPS.rows, FPS.cols, 7],[ 1 -1])
    # grid.f_prop[FPS.rows, FPS.cols, 8] = circshift(grid.f_temp[FPS.rows, FPS.cols, 8],[-1 -1])
    # grid.f_prop[FPS.rows, FPS.cols, 9] = circshift(grid.f_temp[FPS.rows, FPS.cols, 9],[-1  1])

    for i = 1:grid.width 
        i_n = get_next_index(grid, i, grid.width)
        i_p = get_prev_index(grid, i, grid.width)

        for j = 1:grid.length 
            j_n = get_next_index(grid, j, grid.length)
            j_p = get_prev_index(grid, j, grid.length)

            grid.lattices[i, j].f_prop = grid.lattices[i, j].f_temp[1]
            grid.lattices[i, j_p].f_prop = grid.lattices[i, j].f_temp[2]
            grid.lattices[i_p, j].f_prop = grid.lattices[i, j].f_temp[3]
            grid.lattices[i, j_n].f_prop = grid.lattices[i, j].f_temp[4]
            grid.lattices[i_n, j].f_prop = grid.lattices[i, j].f_temp[5]
            grid.lattices[i_p, j_p].f_prop = grid.lattices[i, j].f_temp[6]
            grid.lattices[i_p, i_n].f_prop = grid.lattices[i, j].f_temp[7]
            grid.lattices[i_n, j_n].f_prop = grid.lattices[i, j].f_temp[8]
            grid.lattices[i_n, j_p].f_prop = grid.lattices[i, j].f_temp[9]
        end
    end
end

immutable PressurePeriodicStream_2D{T <: Direction,
                           V <: _2D} <: Streaming

    rho_inlet::Float64
    rho_outlet::Float64
    inlet_row::Int64
    inlet_col::Array{Int64, 1}
    outlet_row::Int64
    outlet_col::Array{Int64, 1}
    indices::Tuple{Base.OnTo{Int64}}

    PressurePeriodicStream_2D(rho_inlet, rho_outlet, 
        inlet_row, inlet_col, outlet_row, outlet_col) = (
            @assert (size(inlet_col) == size(outlet_col));

            new(rho_inlet, rho_outlet, inlet_row, inlet_col,
                outlet_row, outlet_col, indices([inlet_row, inlet_col]))
        )

end

# ============================================================
# ==== Periodic Streaming with Pressure condition
# ============================================================
function periodic_pressure!(grid::Grid_2D, d2q9::D2Q9,
                           bound_in_row::Int64,
                           bound_in_col::Array{Int64, 1},
                           bound_out_row::Int64,
                           bound_out_col::Array{Int64,1},
                           bound_rho::Float64)

    grid.lattices[bound_in_row, bound_in_col].f_temp = 
        f_eq(d2q9, grid.lattices[bound_out_row, bound_out_col], bound_rho)
             (grid.lattices[bound_out_row, bound_out_col].f_temp .-  
             grid.lattices[bound_out_row, bound_out_col].f_eq)

end 

function streaming!(PFPS::PressurePeriodicStream_2D{West, D2Q9}, grid::Grid_2D,
                    d2q9::D2Q9)

    # Compute the densities for the inlet and outlet,
    # where the pressure is known
    #Inlet
    # grid.f_temp[PFPS.inlet_row, PFPS.inlet_col, :] =
        periodic_pressure!(grid, d2q9, 
             PFPS.inlet_row, PFPS.inlet_col ,
             PFPS.outlet_row-1, PFPS.outlet_col ,
                PFPS.rho_inlet)

    #Outlet
    grid.f_temp[PFPS.outlet_row, PFPS.outlet_col, :] =
        periodic_pressure(grid, d2q9, PFPS.inlet_row+1,
                          PFPS.inlet_col,
                          PFPS.rho_outlet)

end

export
    FullPeriodicStreaming_2D,
    PressurePeriodicStream_2D
