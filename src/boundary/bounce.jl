#! /usr/bin/env julia

# ======= Bounces
type Bounce{T <: Direction, V<:Velocity_Set} <: Boundary
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}
end

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
