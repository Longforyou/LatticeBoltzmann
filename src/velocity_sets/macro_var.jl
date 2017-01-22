#! /usr/bin/env julia

"""
This files contains the description of come generic functions
needed for the computation of the macro
"""
function density(f_prop::Array{Float64, 1})
    return sum(f_prop)
end

function compute_macro_var(grid::Grid_2D{_2D, Flow}) 

    for i in 1:grid.width, j in 1:grid.length
         grid.density[i,j] = density( grid.f_prop[i,j, :])
         grid.velocity[i, j, :] = _D2Q9.velo(grid.f_prop[i,j,:])

    end

end

@acc function compute_macro_var(grid::Grid_2D{D2Q9, Compressible})

    for i in 1:grid.width, j in 1:grid.length
         grid.density[i,j] = density( grid.f_prop[i,j, :])
         grid.velocity[i, j, :] = _D2Q9.velo(grid.f_prop[i,j,:])

    end

end

function compute_macro_var(grid::Grid_2D{D2Q9, Incompressible})

    for i in 1:grid.width, j in 1:grid.length
         grid.density[i,j] = density( grid.f_prop[i,j, :])
         grid.velocity[i, j, :] = _D2Q9.velo(grid.f_prop[i,j,:])

    end

end

