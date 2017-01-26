#! /usr/bin/env julia

__precompile__()

module LatticeBoltzmann

using  ProgressMeter, WriteVTK #, ParallelAccelerator


# Process all other files
include("setup.jl")

export compute!

end # module
