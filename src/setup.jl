#!/usr/bin/env julia

# Process the abstract modelling
include("abstract_lbm.jl")
using .Abstract_LBM

include("constants/constants.jl")
using .constants

include("io/output.jl")
include("velocity_sets/velocity_set.jl")
using ._D2Q9.D2Q9

include("streaming/streaming.jl")
include("collision/collision.jl")
include("analytical/analytical.jl")
include("flow/flow.jl")
include("grid/grid.jl")
include("boundary/boundary.jl")

