#!  /usr/bin/env/julia

module _Lattice
using ..Abstract_LBM

type Lattice <: Particle
    
    f_prop::Array{Float64, 1}
    f_eq::Array{Float64, 1}
    f_temp::Array{Float64, 1}
    density::Float64
    velocity::Array{Float64, 1}

    Lattice(directions::Int64, velocities::Int64) = (
        f_prop = zeros(directions);
        f_eq = zeros(directions);
        f_temp = zeros(directions);
        density = 0.;
        velocity = zeros(velocities);

        new(f_prop, f_eq, f_temp, density, velocity);
    )

end 

function get_macro_var(lattice_arr::Array{Lattice})

    sz_latt_arr = size(lattice_arr)
    density_arr = zeros(Float64, sz_latt_arr)

    # TODO: maybe another to generate these matrices 
    velocity_arr_x = zeros(Float64, sz_latt_arr)
    velocity_arr_y = zeros(Float64, sz_latt_arr)

    index = 1;
    for lattice in lattice_arr

        density_arr[index] = lattice.density
        
        # Velocities are swapped! Since julia uses column major formats..
        velocity_arr_x[index] = lattice.velocity[2] / lattice.density
        velocity_arr_y[index] = lattice.velocity[1] / lattice.density

        index += 1
    end

    return density_arr, velocity_arr_x, velocity_arr_y
end

function valid_lattice_direction(direction::Int64)

    curr_direction = direction == 1 ? 2 : (direction == 2 ? 1 : 0)
    assert(curr_direction > 0) 

    return curr_direction
end

function get_lattice_velocity(lattice_arr::Array{Lattice}, row::Array{Int64,1}, col::Array{Int64, 1}, direction::Int64, col_squeeze::Bool=true)

    # Get the correct direction
    curr_direction = valid_lattice_direction(direction)

    out_arr = zeros((length(row), length(col)))

    for i in 1:length(row), j in 1:length(col)
        out_arr[i, j] = lattice_arr[i, j].velocity[direction] ./ lattice_arr[i, j].density
    end

    if col_squeeze
        return squeeze(out_arr, 1)
    else
        return out_arr
    end
end 

export Lattice, get_macro_var, get_lattice_velocity

end
