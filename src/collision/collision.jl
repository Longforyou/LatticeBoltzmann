#! /usr/bin/env julia

"""
Type for the description of the BGK collision operator.
"""
immutable BGK <: Collision
    omega::Float64

    function BGK(consts::LBM_Constants)
        omega = 1. / consts.tau;
        println("Omega: ", omega)
        new(omega)
    end
end

"""
    compute_collision(grid, collision)

Interface function for the computation of different
collision operators via a dispatch function `collide`.
"""
function compute_collision!(grid::Grid, collision::Collision)

    # The collision of the particles is indepentend of
    # the particle type
    grid.f_temp = collide(collision, grid)

end

"""
    collide(bgk, grid)

The BGK-Operator relaxates the population `f_prop`
 to a equilibrium function `f_eq`.
"""
function collide(bgk::BGK, grid::Grid_2D)
    return grid.f_prop .* (1. - bgk.omega) .+
        grid.f_eq .* bgk.omega
end

# Make them available
export BGK
