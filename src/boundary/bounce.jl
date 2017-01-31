#! /usr/bin/env julia

"""
    Bounce{T<:Direction, V<:Velocity} <: Boundary

Wrapper for the indices of a Bounceback condition on the domain.
"""
type Bounce{T <: Direction, V<:Velocity_Set} <: Boundary
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}
end

"""
    boundary(grid, bound::Bounce{T<:Direction, D2Q9})

Implementation of the interface function for the
computation of the bounce boundary condition. Dispatch
over the direction.
"""
function boundary!(grid::Grid_2D,
                  bound::Bounce{North, D2Q9},
                  d2q9::D2Q9)
    
    grid.f_prop[bound.rows, bound.cols, d2q9.dict[South]] =
      grid.f_temp[bound.rows, bound.cols, d2q9.dict[North]]
end

function boundary!(grid::Grid_2D,
                  bound::Bounce{South, D2Q9},
                  d2q9::D2Q9)
    grid.f_prop[bound.rows, bound.cols, d2q9.dict[North]] =
      grid.f_temp[bound.rows, bound.cols, d2q9.dict[South]]
end

function boundary!(grid::Grid_2D,
                  bound::Bounce{West, D2Q9},
                  d2q9::D2Q9)
    grid.f_prop[bound.rows, bound.cols, d2q9.dict[East]] =
      grid.f_temp[bound.rows, bound.cols, d2q9.dict[West]]
end

function boundary!(grid::Grid_2D,
                  bound::Bounce{East, D2Q9},
                  d2q9::D2Q9)
    grid.f_prop[bound.rows, bound.cols, d2q9.dict[West]] =
      grid.f_temp[bound.rows, bound.cols, d2q9.dict[East]]
end
