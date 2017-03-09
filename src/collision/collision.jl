#! /usr/bin/env julia

"""
This module contains the description of collision operators.
Eg. the Bhatnagar-Gross-Kroog operator
"""

immutable BGK <: Collision
    omega::Float64

    function BGK(consts::LBM_Constants)
        omega = 1. / consts.tau;
        println("Omega: ", omega)
        new(omega)
    end
end

# ===============
function compute_collision!(grid::Grid_2D, collision::Collision)

    # The collision of the particles is indepentend of
    # the particle type
    collide!(collision, grid.lattices)

end

# """
# The BGK-Operator relaxates the population 'f_prop'
#  to a equilibrium function.
# """
function collide!(bgk::BGK, lattice::Array{Lattice})

    for index = size(lattice)
        
        lattice[index].f_temp = lattice[index].f_prop .* (1. - bgk.omega) .+ lattice[index].f_eq .* bgk.omega
    end
end

# Make them available
export BGK
