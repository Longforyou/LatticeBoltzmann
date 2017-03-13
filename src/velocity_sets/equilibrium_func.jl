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

    println("Pre EQ:", grid_lattices[3,3].f_eq)
    for i = 1:sz_grid[1]
        for j = 1:sz_grid[2]
            grid_lattices[i, j].f_eq = 
                f_eq_kernel!(grid_lattices[i, j], d2q9)
        end # j
    end # i
    println("Post EQ:", grid_lattices[3,3].f_eq)
end # f_eq

# For the boundary compuations
function f_eq(lattice::Lattice, velset::_2D, rho::Float64)

    f_eq_kernel!(lattice, velset, rho)

    return lattice.f_eq
end # f_eq

function f_eq_kernel!(lattice::Lattice, d2q9::D2Q9{Compressible})

    ueqxij = lattice.velocity[1]; ueqyij = lattice.velocity[2];
    uxsq = ueqxij ^ 2; uysq = ueqyij ^ 2;
    rho = lattice.density

    uxuy6 = ueqxij + ueqyij; uxuy7 = -ueqxij + ueqyij;
    uxuy8 = -ueqxij - ueqyij; uxuy9 = ueqxij - ueqyij;

    usq = 1.5 * (uxsq + uysq)

    lattice.f_eq[1] = copy(d2q9.w[1] * (rho - usq))
    lattice.f_eq[2] = copy(d2q9.w[2] * (rho + 3. * ueqxij + 4.5 * uxsq - usq))
    lattice.f_eq[3] = copy(d2q9.w[3] * (rho + 3. * ueqyij + 4.5 * uysq - usq))
    lattice.f_eq[4] = copy(d2q9.w[4] * (rho - 3. * ueqxij + 4.5 * uxsq - usq))
    lattice.f_eq[5] = copy(d2q9.w[5] * (rho - 3. * ueqyij + 4.5 * uysq - usq))
    lattice.f_eq[6] = copy(d2q9.w[6] * (rho + 3. * uxuy6 + 4.5 * uxuy6^2 - usq))
    lattice.f_eq[7] = copy(d2q9.w[7] * (rho + 3. * uxuy7 + 4.5 * uxuy7^2 - usq))
    lattice.f_eq[8] = copy(d2q9.w[8] * (rho + 3. * uxuy8 + 4.5 * uxuy8^2 - usq))
    lattice.f_eq[9] = copy(d2q9.w[9] * (rho + 3. * uxuy9 + 4.5 * uxuy9^2 - usq))

end

function f_eq_kernel!(lattice::Lattice, d2q9::D2Q9{Incompressible})

    ueqxij = lattice.velocity[1]; ueqyij = lattice.velocity[2];
    uxsq = ueqxij ^ 2; uysq = ueqyij ^ 2;
    rho = lattice.rho

    uxuy6 = ueqxij + ueqyij; uxuy7 = -ueqxij + ueqyij;
    uxuy8 = -ueqxij - ueqyij; uxuy9 = ueqxij - ueqyij;


    lattice.f_eq[1] = copy(d2q9.w[1] * rho)
    lattice.f_eq[2] = copy(d2q9.w[2] * (rho + 3. * ueqxij))
    lattice.f_eq[3] = copy(d2q9.w[3] * (rho + 3. * ueqyij))
    lattice.f_eq[4] = copy(d2q9.w[4] * (rho - 3. * ueqxij))
    lattice.f_eq[5] = copy(d2q9.w[5] * (rho - 3. * ueqyij))
    lattice.f_eq[6] = copy(d2q9.w[6] * (rho + 3. * uxuy6))
    lattice.f_eq[7] = copy(d2q9.w[7] * (rho + 3. * uxuy7))
    lattice.f_eq[8] = copy(d2q9.w[8] * (rho + 3. * uxuy8))
    lattice.f_eq[9] = copy(d2q9.w[9] * (rho + 3. * uxuy9))

end

function f_eq_kernel!(lattice::Lattice, d2q9::D2Q9{Compressible}, rho::Float64)

    ueqxij = lattice.velocity[1]; ueqyij = lattice.velocity[2];
    uxsq = ueqxij ^ 2; uysq = ueqyij ^ 2;

    uxuy6 = ueqxij + ueqyij; uxuy7 = -ueqxij + ueqyij;
    uxuy8 = -ueqxij - ueqyij; uxuy9 = ueqxij - ueqyij;

    usq = 1.5 * (uxsq + uysq)

    lattice.f_eq[1] = copy(d2q9.w[1] * (rho - usq))
    lattice.f_eq[2] = copy(d2q9.w[2] * (rho + 3. * ueqxij + 4.5 * uxsq - usq))
    lattice.f_eq[3] = copy(d2q9.w[3] * (rho + 3. * ueqyij + 4.5 * uysq - usq))
    lattice.f_eq[4] = copy(d2q9.w[4] * (rho - 3. * ueqxij + 4.5 * uxsq - usq))
    lattice.f_eq[5] = copy(d2q9.w[5] * (rho - 3. * ueqyij + 4.5 * uysq - usq))
    lattice.f_eq[6] = d2q9.w[6] * (rho + 3. * uxuy6 + 4.5 * uxuy6^2 - usq)
    lattice.f_eq[7] = d2q9.w[7] * (rho + 3. * uxuy7 + 4.5 * uxuy7^2 - usq)
    lattice.f_eq[8] = d2q9.w[8] * (rho + 3. * uxuy8 + 4.5 * uxuy8^2 - usq)
    lattice.f_eq[9] = d2q9.w[9] * (rho + 3. * uxuy9 + 4.5 * uxuy9^2 - usq)

end

function f_eq_kernel!(lattice::Lattice, d2q9::D2Q9{Incompressible}, rho::Float64)

    ueqxij = lattice.velocity[1]; ueqyij = lattice.velocity[2];
    uxsq = ueqxij ^ 2; uysq = ueqyij ^ 2;

    uxuy6 = ueqxij + ueqyij; uxuy7 = -ueqxij + ueqyij;
    uxuy8 = -ueqxij - ueqyij; uxuy9 = ueqxij - ueqyij;

    lattice.f_eq[1] = d2q9.w[1] * rho
    lattice.f_eq[2] = d2q9.w[2] * (rho + 3. * ueqxij)
    lattice.f_eq[3] = d2q9.w[3] * (rho + 3. * ueqyij)
    lattice.f_eq[4] = d2q9.w[4] * (rho - 3. * ueqxij)
    lattice.f_eq[5] = d2q9.w[5] * (rho - 3. * ueqyij)
    lattice.f_eq[6] = d2q9.w[6] * (rho + 3. * uxuy6)
    lattice.f_eq[7] = d2q9.w[7] * (rho + 3. * uxuy7)
    lattice.f_eq[8] = d2q9.w[8] * (rho + 3. * uxuy8)
    lattice.f_eq[9] = d2q9.w[9] * (rho + 3. * uxuy9)

end
