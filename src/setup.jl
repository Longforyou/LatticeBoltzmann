#!/usr/bin/env julia

"""
This file contains the inclusion of all paths
"""
include("./constants/constants.jl")
include("./boundary/boundary.jl")
include("./velocity_sets/velocity_set.jl")
include("./grid/grid.jl")
include("./io/output.jl")
include("./streaming/streaming.jl")
include("./collision/collision.jl")
include("./analytical/analytical.jl")
include("./flow/flow.jl")
