#! /usr/bin/env julia

"""
This files contains the description of come generic functions
needed for the computation of the macro
"""
function density(f_prop::Array{Float64, 1})
    return sum(f_prop)
end

function compute_macro_var(grid::Grid_2D, d2q9::D2Q9)

    for i = 1:grid.width
        for j = 1:grid.length
         @inbounds @fastmath grid.density[i, j] = density(grid.f_prop[i, j, :])
         @inbounds @fastmath grid.velocity[i, j, 1] = _D2Q9.velo_1(grid.f_prop[i, j, :])
         @inbounds @fastmath grid.velocity[i, j, 2] = _D2Q9.velo_2(grid.f_prop[i, j, :])
        end
    end

end

