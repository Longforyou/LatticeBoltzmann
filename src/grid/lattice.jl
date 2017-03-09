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

export Lattice, get_macro_var

end
