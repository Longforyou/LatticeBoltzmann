#! /usr/bin/env julia

"""
    density(f_prop)

Computes the density for the passed discrete
distribution array of arbitrary length.
"""
function density(f_prop::Array{Float64, 1})
    return sum(f_prop)
end

"""
    compute_macro_var!(grid, velset)

Interface function for the computation of all makroskopic
variables of the grid. Note that the velocity is computed via

\rho\boldsymbol{u} = \sum_i \boldsymbol{c}_i\dot f_i

So the computes velocity contains the density. It's a numerical
optimization of the code (division is expensive).
"""
function compute_macro_var!(grid::Grid_2D, velset::_2D)

    for i = 1:grid.width
        for j = 1:grid.length
         @inbounds @fastmath grid.density[i, j] = density(grid.f_prop[i, j, :])
         @inbounds @fastmath grid.velocity[i, j, 1] = _D2Q9.velo_1(grid.f_prop[i, j, :])
         @inbounds @fastmath grid.velocity[i, j, 2] = _D2Q9.velo_2(grid.f_prop[i, j, :])
        end
    end

end

