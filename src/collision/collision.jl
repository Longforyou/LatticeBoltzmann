#! /usr/bin/env julia

module _Collision

"""
This module contains the description of collision operators.
Eg. the Bhatnagar-Gross-Kroog operator
"""

abstract BGK <: Collision

# ===============
function compute_collision(lbm::LBM{V <: Velocity_Set, F <: Flow,
                                    S <: Streaming, C <: Collision})

    # The collision of the particles is indepentend of
    # the particle type
    lbm.grid.f_temp = collision{C}(lbm.grid.f_prop,
                                   lbm.constants.omega,
                                    lbm.grid.f_eq)

end


"""
The BGK-Operator relaxates the population 'f_prop' to a equilibrium function.
"""
function collision{BGK}(LBM{V, F, S, BGK})
    LBM.grid.f_prop .* (1. - LBM.constants.omega) .+
        LBM.grid.f_eq .* LBM.constants.omega
end

end


