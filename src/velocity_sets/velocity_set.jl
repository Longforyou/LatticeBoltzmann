#! /usr/bin/env julia

# using .Abstract_LBM

"""
This module contains definitions of types for referencing
the different kinds of velocity fields.
Example
-------
using velocity_set._D2Q9

make the weights w and the direction arrays c_x & c_y
 visible in the global namespace
"""

module _D2Q9
  using ..Abstract_LBM
  using .._Lattice

immutable D2Q9{F<:Flow} <: _2D
    
    c_x::Array{Float64, 1}
    c_y::Array{Float64, 1}
    w::Array{Float64, 1}
    dict::Dict{DataType, Array{Int64, 1}}

    D2Q9() = (
      # Some constant fields
      new(Array{Float64, 1}([0., 1., 0., -1., 0., 1., -1., -1., 1.]), 
      Array{Float64, 1}([0., 0., 1., 0., -1., 1., 1., -1., -1.]),
      Array{Float64, 1}([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]),
          Dict{DataType, Array{Int64, 1}}([(North, [3, 6, 7]), (South, [5, 8, 9]),
                (West, [4, 7, 8]), (East, [2, 6, 9])]))
    )

  end

function velo!(lattice::Lattice)
  lattice.velocity[1] = velo_1(lattice.f_prop)  
  lattice.velocity[2] = velo_2(lattice.f_prop) 
  _Lattice.set_velocity(lattice, [velo_1(lattice.f_prop), velo_2(lattice.f_prop)])
end

function velo_1(f_prop::Array{Float64, 1})

    return f_prop[2] + f_prop[6] + f_prop[9] -
        f_prop[4] - f_prop[7] - f_prop[8]

end

function velo_2(f_prop::Array{Float64, 1})
    return f_prop[3] + f_prop[6] + f_prop[7] -
        f_prop[5] - f_prop[8] - f_prop[9]
end

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

  # Make the variables visible in the global namespace
  export D2Q9

end # module _D2Q9

# Include the 2D velocity sets
# include("velocity_2d.jl"); export _D2Q9

# TODO 3d velocity sets

using ._D2Q9

include("equilibrium_func.jl")
include("macro_var.jl")

