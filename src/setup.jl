#!/usr/bin/env julia

# Process the abstract modelling
include("abstract_lbm.jl")
using .Abstract_LBM

include("./constants/constants.jl")
include("./io/output.jl")
include("./streaming/streaming.jl")
include("./collision/collision.jl")
include("./analytical/analytical.jl")
include("./flow/flow.jl")
include("./velocity_sets/velocity_set.jl")
include("./grid/grid.jl")
include("./boundary/boundary.jl")

using velocity_set, velocity_set._D2Q9
