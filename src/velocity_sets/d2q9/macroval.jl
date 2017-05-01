#! usr/bin/julia

function velo!(lattice::Lattice, d2q9::D2Q9)
    lattice.velocity[1] = firstVelocity(lattice.f_prop, d2q9)
    lattice.velocity[2] = secondVelocity(lattice.f_prop, d2q9)
end

"""
    firstVelocity(f_prop, d2q9)

Computes the first velocity component for the passed
9 valued distribution function.
"""
function firstVelocity(f_prop::Array{Float64, 1}, d2q9::D2Q9)

    return f_prop[2] + f_prop[6] + f_prop[9] -
        f_prop[4] - f_prop[7] - f_prop[8]

end

"""
    secondVelocity(f_prop, d2q9)

Computes the second velocity component for the passed
9 valued distribution function.
"""
function secondVelocity(f_prop::Array{Float64, 1}, d2q9::D2Q9)
    return f_prop[3] + f_prop[6] + f_prop[7] -
        f_prop[5] - f_prop[8] - f_prop[9]
end
