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
    # println("Collision")
    # print_lattice_f_temp(grid.lattices)
    collide!(collision, grid.lattices)

    # println("Post")
    # print_lattice_f_temp(grid.lattices)

end

"""
    collide(bgk, grid)

The BGK-Operator relaxates the population `f_prop`
 to a equilibrium function `f_eq`.
"""
function collide!(bgk::BGK, lattice::Array{Lattice})

    for index = 1:length(lattice)
        _Lattice.set_f_temp!(lattice[index], lattice[index].f_prop .* (1. - bgk.omega) .+ lattice[index].f_eq .* bgk.omega)
    end
end

# Make them available
export BGK
