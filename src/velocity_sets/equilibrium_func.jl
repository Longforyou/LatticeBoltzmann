#! /usr/bin/env julia

"""
    compute_f_eq!(grid, velset)

Compute the equilibrium equation for each node in the
grid. 
"""
function compute_f_eq!(grid::Grid_2D, velset::_2D)
    f_eq!(grid.lattices, velset)
end

include("d2q9/equilibrium.jl")

