#! /usr/bin/env julia

"""
This module contains the description of collision operators.
Eg. the Bhatnagar-Gross-Kroog operator
"""

abstract BGK <: Collision

# ===============
function compute_collision(lbm::LBM{Velocity_Set, Flow,
                                    Streaming, Collision})

    # The collision of the particles is indepentend of
    # the particle type
    lbm.grid.f_temp = collision(lbm)

end


"""
The BGK-Operator relaxates the population 'f_prop' to a equilibrium function.
"""
function collision(lbm::LBM{Velocity_Set, Flow,
                       Streaming, BGK})
    lbm.grid.f_prop .* (1. - lbm.constants.omega) .+
        lbm.grid.f_eq .* lbm.constants.omega
end



