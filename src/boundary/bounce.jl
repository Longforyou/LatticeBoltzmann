#! /usr/bin/env julia

# ======= Bounces
type Bounce{T <: Direction, V<:Velocity_Set} <: Boundary
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}
end

function boundary(grid::Grid_2D,
                  bound::Bounce{North, D2Q9})
    
    grid.f_prop[bound.rows, bound.cols, [4, 7, 8]] =
      grid.f_prop[bound.rows, bound.cols, [2, 5, 6]]
end

function boundary(grid::Grid_2D,
                  bound::Bounce{South, D2Q9})
    grid.f_prop[bound.rows, bound.cols, [2, 5, 6]] =
      grid.f_prop[bound.rows, bound.cols, [4, 7, 8]]
end

function boundary(grid::Grid_2D,
                  bound::Bounce{West, D2Q9})
    grid.f_prop[bound.rows, bound.cols, [3, 7, 6]] =
      grid.f_prop[bound.rows, bound.cols, [1, 5, 8]]
end

function boundary(grid::Grid_2D,
                  bound::Bounce{East, D2Q9})
    grid.f_prop[bound.rows, bound.cols, [1, 5, 8]] =
      grid.f_prop[bound.rows, bound.cols, [3, 7, 6]]
end
