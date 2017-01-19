#! /usr/bin/env julia

# Open Boundary
type OpenBound{T <: Direction, V<:Velocity_Set._2D} <: Boundary
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}
 end

# ============================================================
# === Open Boundary
# ============================================================
function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::OpenBound{East})
    lbm.grid.f_prop[bound.rows, bound.cols, [1, 5, 8]] =
        2. * lbm.grid.f_prop[bound.rows-1, bound.cols, [1, 5, 8]] -
        lbm.grid.f_prop[bound.rows-2, bound.cols, [1, 5, 8]]

    lbm.grid.velocity[bound.rows, bound.cols, 2] = 0.
end
