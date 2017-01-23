#! /usr/bin/env julia

"""
This module contains the description of collision operators.
Eg. the Bhatnagar-Gross-Kroog operator
"""

immutable BGK <: Collision
    omega::Float64

    function BGK(consts::LBM_Constants)
        omega = 1. / consts.tau
        new(omega)
    end
end

# ===============
function compute_collision(grid::Grid_2D, collision::BGK)

    # The collision of the particles is indepentend of
    # the particle type
    grid.f_temp = collision(grid)

end

"""
The BGK-Operator relaxates the population 'f_prop'
 to a equilibrium function.
"""
function (bgk::BGK)(grid::Grid_2D)
    return grid.f_prop .* (1. - bgk.omega) .+
        grid.f_eq .* bgk.omega
end


# Make them available
export BGK
