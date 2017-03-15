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

    # println("Streaming")
    # Distribution direction
    # println("Pre \n")
    # print_lattice_f_prop(grid.lattices)

    for i = 1:grid.width 
        # println("i_p \t i \t i_n\n")
        i_n = get_next_index(grid, i, grid.width)
        i_p = get_prev_index(grid, i, grid.width)
        # println(i_p, "\t", i, "\t", i_n)

        for j = 1:grid.length 
            # println("j_p \t j \t j_n\n")
            j_n = get_next_index(grid, j, grid.length)
            j_p = get_prev_index(grid, j, grid.length)
            
            #println(j_p, "\t", j, "\t", j_n)

            _Lattice.set_f_prop!(grid.lattices[i, j], 1, grid.lattices[i, j].f_temp[1])
            _Lattice.set_f_prop!(grid.lattices[i_n, j], 2, grid.lattices[i, j].f_temp[2])
            _Lattice.set_f_prop!(grid.lattices[i, j_n], 3, grid.lattices[i, j].f_temp[3])
            _Lattice.set_f_prop!(grid.lattices[i_p, j], 4, grid.lattices[i, j].f_temp[4])
            _Lattice.set_f_prop!(grid.lattices[i, j_p], 5, grid.lattices[i, j].f_temp[5])
            _Lattice.set_f_prop!(grid.lattices[i_n, j_n], 6, grid.lattices[i, j].f_temp[6])
            _Lattice.set_f_prop!(grid.lattices[i_p, j_n], 7, grid.lattices[i, j].f_temp[7])
            _Lattice.set_f_prop!(grid.lattices[i_p, j_p], 8, grid.lattices[i, j].f_temp[8])
            _Lattice.set_f_prop!(grid.lattices[i_n, j_p], 9, grid.lattices[i, j].f_temp[9])
        end
    end
    # println("Post \n")
    # print_lattice_f_prop(grid.lattices)
end

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
