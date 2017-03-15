#! /usr/bin/env julia

"""
This files contains the description of come generic functions
needed for the computation of the macro
"""
function density(f_prop::Array{Float64, 1})
    return sum(f_prop)
end

function compute_macro_var!(grid::Grid_2D, d2q9::D2Q9)

    for i = 1:grid.width, j = 1:grid.length
         @inbounds @fastmath _Lattice.set_density!(grid.lattices[i, j],
            density(grid.lattices[i, j].f_prop))
         @inbounds @fastmath  _D2Q9.velo!(grid.lattices[i, j])
            print_lattice_macro_var(grid.lattices[i, j])
    end



end

