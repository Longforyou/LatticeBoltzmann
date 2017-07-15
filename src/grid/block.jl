#! usr/bin/julia

"""
    Block{D <: UInt8, Q <: UInt8}

A block contains the two arrays associated with every LBM simulation,
the macroscopic variables and the lattice populations.
"""
mutable struct Block{D <: UInt8, Q <: UInt8} <: AbstractBlock
    macroVar::Array{Float64, D}
    populations::Array{Float64, Q}

end

function computeNewWeight(block::Block, maxLvL::Int64)
    return 2 ^ (block.weight - maxLvL)

end
