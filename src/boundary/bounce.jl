#! /usr/bin/env julia

# ======= Bounces
type Bounce{T <: Direction, V<:Velocity_Set} <: Boundary
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}
end

function boundary!(grid::Grid_2D,
                  bound::Bounce{North, D2Q9},
                  d2q9::D2Q9)
    
    grid.lattice[bound.rows, bound.cols].lattice[d2q9dict[South]] =
      grid.lattice[bound.rows, bound.cols].f_temp[d2q9dict[North]]
end

function boundary!(grid::Grid_2D,
                  bound::Bounce{South, D2Q9},
                  d2q9::D2Q9)
    grid.lattice[bound.rows, bound.cols].lattice[d2q9dict[North]] =
      grid.lattice[bound.rows, bound.cols].f_temp[d2q9dict[South]]
end

function boundary!(grid::Grid_2D,
                  bound::Bounce{West, D2Q9},
                  d2q9::D2Q9)
    grid.lattice[bound.rows, bound.cols].f_prop[dict[East]] =
      grid.lattice[bound.rows, bound.cols].f_temp[dict[West]]
end

function boundary!(grid::Grid_2D,
                  bound::Bounce{East, D2Q9},
                  d2q9::D2Q9)
    grid.lattice[bound.rows, bound.cols].f_prop[dict[West]] =
      grid.lattice[bound.rows, bound.cols].f_temp[dict[East]]
end
