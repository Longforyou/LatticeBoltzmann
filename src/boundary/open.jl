#! /usr/bin/env julia

# Open Boundary
type OpenBound{T <: Direction, V<:Velocity_Set} <: Boundary
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}
end

# ============================================================
# === Open Boundary
# ============================================================
function boundary(grid::Grid_2D,
                  bound::OpenBound{East, D2Q9})

    grid.f_prop[bound.rows, bound.cols, [1, 5, 8]] =
        2. * grid.f_prop[bound.rows-1, bound.cols, [1, 5, 8]] -
        grid.f_prop[bound.rows-2, bound.cols, [1, 5, 8]]

    grid.velocity[bound.rows, bound.cols, 2] = 0.
end
