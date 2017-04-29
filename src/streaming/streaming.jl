#! /usr/bin/env julia

"""
    FullPeriodicStreaming_2D

Streaming operator for the streaming of the entire domain.
"""
immutable FullPeriodicStreaming_2D <: Streaming
    rows::UnitRange{Int64}
    cols::UnitRange{Int64}

    function FullPeriodicStreaming_2D(grid::Grid_2D)
        new(UnitRange{Int64}(1:grid.width),
            UnitRange{Int64}(1:grid.length))
    end

end

# =========== Streaming
abstract InnerStreaming <: Streaming #TODO implement different operators..

"""
    compute_streaming!(grid, stream, velset)

Interface function for the computation of all passed Array of
Streaming-Operators `stream`. The dispatch is done via the
function `streaming!`.
"""
function compute_streaming!(grid::Grid, stream::Array{Streaming, 1}, velset::Velocity_Set)
    for stre in stream
        @inbounds streaming!(stre, grid, velset)
    end
end

"""
    streaming!(FPS, grid, d2q9)

FullPeriodicStreaming_2D operator definition. Performs circshifts
on different dimensions of the collided particles to the neighbouring
particles.
"""
function streaming!(FPS::FullPeriodicStreaming_2D, grid::Grid_2D,
                    d2q9::D2Q9)

    for i = 1:grid.width 
        i_n = get_next_index(grid, i, grid.width)
        i_p = get_prev_index(grid, i, grid.width)

        for j = 1:grid.length 
            j_n = get_next_index(grid, j, grid.length)
            j_p = get_prev_index(grid, j, grid.length)
            
            grid.lattices[i, j].f_prop[1] =
                grid.lattices[i, j].f_temp[1]
            grid.lattices[i_n, j].f_prop[2] =
                grid.lattices[i, j].f_temp[2]
            grid.lattices[i, j_n].f_prop[3] =
                grid.lattices[i, j].f_temp[3]
            grid.lattices[i_p, j].f_prop[4] =
                grid.lattices[i, j].f_temp[4]
            grid.lattices[i, j_p].f_prop[5] =
                grid.lattices[i, j].f_temp[5]
            grid.lattices[i_n, j_n].f_prop[6] =
                grid.lattices[i, j].f_temp[6]
            grid.lattices[i_p, j_n].f_prop[7] =
                grid.lattices[i, j].f_temp[7]
            grid.lattices[i_p, j_p].f_prop[8] =
                grid.lattices[i, j].f_temp[8]
            grid.lattices[i_n, j_p].f_prop =
                grid.lattices[i, j].f_temp[9]
        end
    end
end

# ============================================================
# ==== Periodic STreaming with Pressure condition
# ============================================================

"""
    PressurePeriodicStream_2D

Contains a description of a periodic pressure/ density gradient.
"""
immutable PressurePeriodicStream_2D{T <: Direction,
                           V <: _2D} <: Streaming

    rho_inlet::Float64
    rho_outlet::Float64
    inlet_row::Int64
    inlet_col::Array{Int64, 1}
    outlet_row::Int64
    outlet_col::Array{Int64, 1}
    length_col::Int64

    PressurePeriodicStream_2D(rho_inlet, rho_outlet, 
        inlet_row, inlet_col, outlet_row, outlet_col) = (
            @assert (size(inlet_col) == size(outlet_col));

            new(rho_inlet, rho_outlet, inlet_row, inlet_col,
                outlet_row, outlet_col, length(inlet_col))
        )

end

# ============================================================
# ==== Periodic Streaming with Pressure condition
# ============================================================
"""
    periodic_pressure(grid, d2q9, bound_row, bound_col, bound_rho)

Function to interface the function `f_eq` in order to compute the value of the discrete distribution function. **NOTE** that the object of the type `PressurePeriodicStream_2D` has to be computed --before-- the streaming of the domain.
"""
function periodic_pressure(grid::Grid_2D, d2q9::D2Q9,
                           bound_in_row::Int64,
                           bound_in_col::Int64,
                           bound_out_row::Int64,
                           bound_out_col::Int64,
                           bound_rho::Float64, directions::Array{Int64, 1})

    
    
    _Lattice.set_f_temp!(grid.lattices[bound_in_row, bound_in_col], 
        f_eq(grid.lattices[bound_out_row, bound_out_col], d2q9, bound_rho) .+
             (grid.lattices[bound_out_row, bound_out_col].f_temp .-  
             grid.lattices[bound_out_row, bound_out_col].f_eq), directions)

end 

"""
    streaming(PFPS, grid, d2q9)

Streaming-Operator implementation for the type `PressurePeriodicStream`.
TODO implement different direction of different directions.
"""
function streaming!(PFPS::PressurePeriodicStream_2D{West, D2Q9}, grid::Grid_2D,
                    d2q9::D2Q9)

    # Compute the densities for the inlet and outlet,
    # where the pressure is known

    # println("Streaming Peri Pres")
    # println("PPre \n")
    # print_lattice_f_temp(grid.lattices[PFPS.inlet_row, PFPS.inlet_col])
    
    for i = 1:PFPS.length_col
        #Inlet
        periodic_pressure(grid, d2q9, 
             PFPS.inlet_row, PFPS.inlet_col[i],
             PFPS.outlet_row-1, PFPS.outlet_col[i] ,
                PFPS.rho_inlet, [1, 5, 8])

        #Outlet
        periodic_pressure(grid, d2q9, 
             PFPS.outlet_row, PFPS.outlet_col[i] ,
             PFPS.inlet_row+1, PFPS.inlet_col[i],
                PFPS.rho_inlet, [3, 6, 7])
    end
    # println("Post \n")
    # print_lattice_f_temp(grid.lattices[PFPS.inlet_row, PFPS.inlet_col])

end

export
    FullPeriodicStreaming_2D,
    PressurePeriodicStream_2D
