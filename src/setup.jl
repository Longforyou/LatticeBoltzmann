#!/usr/bin/env julia

using StaticArrays

# Process the abstract modelling
include("abstract_lbm.jl")
include("constants/constants.jl")
include("geometry/geometry.jl")
include("grid/lattice.jl")
include("grid/grid.jl")
include("io/output.jl")
include("velocity_sets/velocity_set.jl")
include("boundary/boundary.jl")
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
    NonEqBounce,
    OpenBounce,
    PeridicPressure
    North, South,
    West, East, Incompressible,
    Compressible,

    # Functions
    get_pressure_pois_1,
    get_pressure_pois_2,
    get_velo_pois_1,
    get_velo_pois_2,
    get_next_index,
    get_axis_vec



