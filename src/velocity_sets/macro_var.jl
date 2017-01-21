#! /usr/bin/env julia

"""
This files contains the description of come generic functions
needed for the computation of the macro
"""
function density(f_prop::Array{Float64, 1})
    return sum(f_prop)
end

function compute_macro_var(lbm::LBM{_2D, Flow,
                                    Streaming, Collision}) 

    for i in 1:lbm.grid.width, j in 1:lbm.grid.length
        lbm.grid.density[i,j] = density(lbm.grid.f_prop[i,j, :])
        lbm.grid.velocity[i, j, :] = velo(_2D, lbm.grid.f_prop[i,j,:])

    end

end

