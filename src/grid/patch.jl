#! usr/bin/julia

mutable struct Patch <: AbstractPatch
    blocks::Array{Block, 1}
    boundary::Array{Boundary, 1}
    level::Int64 # Level of the patch
    weight::Int64 # Computational weight
    neighbour::Array{Int64, 1} # The indices of the blocks
    adjacent::Array{Int64, 1} # The indices of adjacent patches

end
