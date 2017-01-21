#!/usr/bin/env julia

# Process the abstract modelling
include("abstract_lbm.jl")
importall .Abstract_LBM

include("constants/constants.jl")

include("io/output.jl")
include("grid/grid.jl")
include("velocity_sets/velocity_set.jl")
using ._D2Q9
include("boundary/boundary.jl")

include("flow/flow.jl")
include("streaming/streaming.jl")
include("collision/collision.jl")
include("analytical/analytical.jl")

# Export all values

export
    # Types
    Grid,
    LBM_Constants,
    Collision,
    Neumann,
    Bounce,
    OpenBounce,
    PeridicPressure
    North, South,
    West, East, Incompressible,
    Compressible,
    LBM_Incompressible,

    # Functions
    get_pressure_pois,
    get_velo_pois,
    getNeighbours,
    get_next_index,
    get_axis_vec



