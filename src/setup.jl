#!/usr/bin/env julia

# Process the abstract modelling
include("abstract_lbm.jl")
importall .Abstract_LBM

include("constants/constants.jl")
using .constants

include("io/output.jl")
include("grid/grid.jl")
include("velocity_sets/velocity_set.jl")
using ._D2Q9
include("boundary/boundary.jl")

include("streaming/streaming.jl")
include("collision/collision.jl")
include("analytical/analytical.jl")
include("flow/flow.jl")

