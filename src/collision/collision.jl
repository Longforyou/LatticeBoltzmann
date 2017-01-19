#! /usr/bin/env julia

module _Collision

"""
This module contains the description of collision operators.
Eg. the Bhatnagar-Gross-Kroog operator
"""

abstract BGK <: Collision

"""
The BGK-Operator relaxates the population 'f_prop' to a equilibrium function.
"""
function collision{BGK}(LBM{V, F, S, BGK})
    LBM.grid.f_prop .* (1. - LBM.constants.omega) .+
        LBM.grid.f_eq .* LBM.constants.omega
end

end


