#! usr/bin/julia

"""
    f_eq!(grid_f_eq, velset, rho, velo)

Interface function to the kernel function `f_eq_kernel`.
Performace the actual iteration over the grid.
"""
function f_eq!(grid_lattices::Array{Lattice, 2}, d2q9::D2Q9)

    sz_grid = size(grid_lattices)

    for i = 1:sz_grid[1]
        for j = 1:sz_grid[2]
            _D2Q9.f_eq_kernel!(grid_lattices[i, j], d2q9)

        end # j
    end # i
end # f_eq

"""
    f_eq(d2q9, rho, velo)

Interface function to the kernel function `f_eq_kernel` for
a subset of particles. Useful for computing boudary conditions.
"""
# For the boundary compuations
function f_eq(lattice::Lattice, d2q9::D2Q9, rho::Float64)

    # println("Called f_eq_for bounds")
    eq = zeros(9)
    eq = _D2Q9.f_eq_kernel(lattice.velocity, eq, d2q9, rho)
    # println("Post: ", lattice.f_eq, eq)

    return eq
end # f_eq

function f_eq_kernel!(lattice::Lattice, d2q9::D2Q9)

    _Lattice.set_f_eq!(lattice, f_eq_kernel(lattice.velocity, lattice.f_eq, d2q9, lattice.density))
end

function f_eq_kernel(velocity::Array{Float64, 1}, f_eq::Array{Float64, 1},
    d2q9::D2Q9{Compressible}, rho::Float64)

    #println("Rho: ", rho)
    #println("f_eq:", f_eq)
    #println("lattice.f_temp - lattice.f_eq", lattice.f_temp - lattice.f_eq)

    ueqxij = velocity[1]; ueqyij = velocity[2];
    uxsq = ueqxij ^ 2; uysq = ueqyij ^ 2;

    uxuy6 = ueqxij + ueqyij; uxuy7 = -ueqxij + ueqyij;
    uxuy8 = -ueqxij - ueqyij; uxuy9 = ueqxij - ueqyij;

    usq = 1.5 * (uxsq + uysq)

    #println(d2q9.w[1] * (rho - usq))
    f_eq[1] = copy(d2q9.w[1] * (rho - usq))
    f_eq[2] = copy(d2q9.w[2] * (rho + 3. * ueqxij + 4.5 * uxsq - usq))
    f_eq[3] = copy(d2q9.w[3] * (rho + 3. * ueqyij + 4.5 * uysq - usq))
    f_eq[4] = copy(d2q9.w[4] * (rho - 3. * ueqxij + 4.5 * uxsq - usq))
    f_eq[5] = copy(d2q9.w[5] * (rho - 3. * ueqyij + 4.5 * uysq - usq))
    f_eq[6] = copy(d2q9.w[6] * (rho + 3. * uxuy6 + 4.5 * uxuy6^2 - usq))
    f_eq[7] = copy(d2q9.w[7] * (rho + 3. * uxuy7 + 4.5 * uxuy7^2 - usq))
    f_eq[8] = copy(d2q9.w[8] * (rho + 3. * uxuy8 + 4.5 * uxuy8^2 - usq))
    f_eq[9] = copy(d2q9.w[9] * (rho + 3. * uxuy9 + 4.5 * uxuy9^2 - usq))

    return f_eq

end

function f_eq_kernel!(velocity::Array{Float64, 1}, f_eq::Array{Float64, 1},
  d2q9::D2Q9{Incompressible}, rho::Float64)

    ueqxij = velocity[1]; ueqyij = velocity[2];
    uxsq = ueqxij ^ 2; uysq = ueqyij ^ 2;

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
