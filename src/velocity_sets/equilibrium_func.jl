#! /usr/bin/env julia

# """
# This file contains a description of usefull function
#  for the computation of the equilibrium function
# """
# =========== Equilibrium Distribution function
function compute_f_eq!(grid::Grid_2D, velset::_2D)
    f_eq!(grid.lattices, velset)
end

function f_eq!(grid_lattices::Array{Lattice, 2}, d2q9::D2Q9)

    sz_grid = size(grid_lattices)

    for i = 1:sz_grid[1]
        for j = 1:sz_grid[2]
            _D2Q9.f_eq_kernel!(grid_lattices[i, j], d2q9)

        end # j
    end # i
end # f_eq

# For the boundary compuations
function f_eq(lattice::Lattice, d2q9::D2Q9, rho::Float64)

    # println("Called f_eq_for bounds")
    eq = zeros(9)
    eq = _D2Q9.f_eq_kernel(lattice.velocity, eq, d2q9, rho)
    # println("Post: ", lattice.f_eq, eq)

    return eq
end # f_eq


