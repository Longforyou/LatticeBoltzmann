#! usr/bin/julia

immutable MultiBlock <: AbstractBlock
    blocks::DistributedDict{Int64, Block{UInt8}}


end

immutable Block{N <: UInt8} <: AbstractBlock
    grid::Array{Lattice, N}
    neighbour::Array{Int64}
    level::Int64
    weight::Int64

end
