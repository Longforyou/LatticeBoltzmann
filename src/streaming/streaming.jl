#! /usr/bin/env julia

immutable FullPeriodicStreaming_2D <: Streaming
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}

    function FullPeriodicStreaming_2D(grid::Grid_2D)
        new(Array{Int64, 1}(1:grid.width),
            Array{Int64, 1}(1:grid.length))
    end

end

abstract InnerStreaming <: Streaming

# =========== Streaming
function compute_streaming(grid::Grid_2D, stream::Array{Streaming, 1})
    for stre in stream
        stre(grid)
    end
end

function (FPS::FullPeriodicStreaming_2D)(grid::Grid_2D)

    # Distribution direction
    grid.f_prop[FPS.rows, FPS.cols, 1] =
        circshift( grid.f_temp[:, : ,1],[ 0  1])
    grid.f_prop[FPS.rows, FPS.cols, 2] =
        circshift( grid.f_temp[FPS.rows, FPS.cols, 2],[ 1  0])
    grid.f_prop[FPS.rows, FPS.cols, 3] =
        circshift( grid.f_temp[FPS.rows, FPS.cols, 3],[ 0 -1])
    grid.f_prop[FPS.rows, FPS.cols, 4] =
        circshift( grid.f_temp[FPS.rows, FPS.cols, 4],[-1  0])
    grid.f_prop[FPS.rows, FPS.cols, 5] =
        circshift( grid.f_temp[FPS.rows, FPS.cols, 5],[ 1  1])
    grid.f_prop[FPS.rows, FPS.cols, 6] =
        circshift( grid.f_temp[FPS.rows, FPS.cols, 6],[ 1 -1])
    grid.f_prop[FPS.rows, FPS.cols, 7] =
        circshift( grid.f_temp[FPS.rows, FPS.cols, 7],[-1 -1])
    grid.f_prop[FPS.rows, FPS.cols, 8] =
        circshift( grid.f_temp[FPS.rows, FPS.cols, 8],[-1  1])
    grid.f_prop[FPS.rows, FPS.cols, 9] =
        grid.f_temp[FPS.rows, FPS.cols, 9]
    
end


#! /usr/bin/env julia

immutable PressurePeriodicStream_2D{T <: Direction,
                           V <: _2D} <: Streaming

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

function (PFPS::PressurePeriodicStream_2D{West, D2Q9})(grid::Grid_2D)

    # Compute the densities for the inlet and outlet,
    # where the pressure is known
    for k in 1:9
        #Inlet
        grid.f_temp[PFPS.inlet_row, PFPS.inlet_col, k] =
            periodic_pressure(grid, k, PFPS.outlet_row-1,
                              PFPS.outlet_col,
                              PFPS.rho_inlet, _D2Q9.w[k], _D2Q9.c_x[k], _D2Q9.c_y[k])

        #Outlet
        grid.f_temp[PFPS.outlet_row, PFPS.outlet_col, k] =
            periodic_pressure(grid, k, PFPS.inlet_row+1,
                              PFPS.inlet_col,
                              PFPS.rho_outlet, _D2Q9.w[k], _D2Q9.c_x[k], _D2Q9.c_y[k])

    end

end

export
    FullPeriodicStreaming_2D,
    PressurePeriodicStream_2D
