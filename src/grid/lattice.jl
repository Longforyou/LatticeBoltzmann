#!  /usr/bin/env/julia

"""
    _Lattice

Contains a description of the type *Lattice*, which is a container for the populations of a lattice boltzmann modell.
All 
"""

immutable Lattice <: Particle
    
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

function set_f_eq!(lattice::Lattice, f_eq::Array{Float64, 1})
    lattice.f_eq = f_eq
end

function set_f_temp!(lattice::Lattice, f_temp::Array{Float64, 1})
    # println("f_temp: \n", f_temp)
    lattice.f_temp = f_temp
end

function set_f_temp!(lattice::Lattice, f_temp::Array{Float64, 1}, directions::Array{Int64, 1})
    # println("f_temp: \n", f_temp)
    lattice.f_temp[directions] = f_temp[directions]
end

function set_f_prop!(lattice::Lattice, i::Int64, f_temp::Float64)
    lattice.f_prop[i] = f_temp
end

function set_density!(lattice::Lattice, d::Float64)
    lattice.density = d
end

function set_velocity!(lattice::Lattice, velo::Array{Float64, 1})
    lattice.velocity = velo
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

function print_lattice_f_eq(lattice_arr::Array{Lattice})
    for i = 1:length(lattice_arr)
        println(lattice_arr[i].f_eq)
    end
    
end 
function print_lattice_f_prop(lattice_arr::Array{Lattice})
    for i = 1:length(lattice_arr)
        println(lattice_arr[i].f_prop)
    end
    
end 
function print_lattice_f_temp(lattice_arr::Array{Lattice})
    for i = 1:length(lattice_arr)
        println(lattice_arr[i].f_temp)
    end
    
end 

function print_lattice_macro_var(lattice::Lattice)
    println("Velocity: ", lattice.velocity ./ lattice.density)
    println("Density: ", lattice.density)
end

export 
    Lattice,
    get_macro_var,
    get_lattice_velocity,
    print_lattice_f_eq,
    print_lattice_f_prop,
    print_lattice_f_temp,
    print_lattice_macro_var
