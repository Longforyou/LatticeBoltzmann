#! /usr/bin/env julia

# """
# This file contains a description of usefull function
#  for the computation of the equilibrium function
# """
# =========== Equilibrium Distribution function
function compute_f_eq(grid::Grid_2D, velset::_2D)
    @inbounds @fastmath f_eq!(grid.f_eq, velset, grid.density, grid.velocity)
end

function f_eq!(grid_f_eq::Array{Float64, 3}, d2q9::D2Q9,
              rho::Array{Float64, 2}, velo::Array{Float64, 3})

    sz_grid = size(grid_f_eq)

    for i = 1:sz_grid[1]
        for j = 1:sz_grid[2]
            f_eq_kernel!(grid_f_eq[i, j, :], rho[i, j], velo[i, j, 1],
                         velo[i, j, 2], d2q9)
        end # j
    end # i
end # f_eq

# For the boundary compuations
function f_eq(d2q9::D2Q9, rho::Float64, velo::Array{Float64, 2})

    sz_grid = size(velo)
    out_f_eq = zeros(sz_grid[1], 9)

    for i = 1:sz_grid[1]
            f_eq_kernel!(out_f_eq[i, :], rho, velo[i, 1], velo[i, 2], d2q9)
    end # i

    return out_f_eq
end # f_eq

function f_eq_kernel!(out_f_eq::Array{Float64, 1}, rho::Float64,
                      ux::Float64, uy::Float64, d2q9::D2Q9{Compressible})

    ueqxij = ux; ueqyij = uy;
    uxsq = ueqxij ^ 2; uysq = ueqyij ^ 2;

    uxuy6 = ueqxij + ueqyij; uxuy7 = -ueqxij + ueqyij;
    uxuy8 = -ueqxij - ueqyij; uxuy9 = ueqxij - ueqyij;

    usq = 1.5 * (uxsq + uysq)

    out_f_eq[1] = d2q9.w[1] * (rho - usq)
    out_f_eq[2] = d2q9.w[2] * (rho + 3. * ueqxij + 4.5 * uxsq - usq)
    out_f_eq[3] = d2q9.w[3] * (rho + 3. * ueqyij + 4.5 * uysq - usq)
    out_f_eq[4] = d2q9.w[4] * (rho - 3. * ueqxij + 4.5 * uxsq - usq)
    out_f_eq[5] = d2q9.w[5] * (rho - 3. * ueqyij + 4.5 * uysq - usq)
    out_f_eq[6] = d2q9.w[6] * (rho + 3. * uxuy6 + 4.5 * uxuy6^2 - usq)
    out_f_eq[7] = d2q9.w[7] * (rho + 3. * uxuy7 + 4.5 * uxuy7^2 - usq)
    out_f_eq[8] = d2q9.w[8] * (rho + 3. * uxuy8 + 4.5 * uxuy8^2 - usq)
    out_f_eq[9] = d2q9.w[9] * (rho + 3. * uxuy9 + 4.5 * uxuy9^2 - usq)

end

function f_eq_kernel!(out_f_eq::Array{Float64, 1}, rho::Float64,
                      ux::Float64, uy::Float64, d2q9::D2Q9{Incompressible})

    ueqxij = ux; ueqyij = uy;
    uxsq = ueqxij ^ 2; uysq = ueqyij ^ 2;

    uxuy6 = ueqxij + ueqyij; uxuy7 = -ueqxij + ueqyij;
    uxuy8 = -ueqxij - ueqyij; uxuy9 = ueqxij - ueqyij;

    out_f_eq[1] = d2q9.w[1] * rho
    out_f_eq[2] = d2q9.w[2] * (rho + 3. * ueqxij)
    out_f_eq[3] = d2q9.w[3] * (rho + 3. * ueqyij)
    out_f_eq[4] = d2q9.w[4] * (rho - 3. * ueqxij)
    out_f_eq[5] = d2q9.w[5] * (rho - 3. * ueqyij)
    out_f_eq[6] = d2q9.w[6] * (rho + 3. * uxuy6)
    out_f_eq[7] = d2q9.w[7] * (rho + 3. * uxuy7)
    out_f_eq[8] = d2q9.w[8] * (rho + 3. * uxuy8)
    out_f_eq[9] = d2q9.w[9] * (rho + 3. * uxuy9)

end
