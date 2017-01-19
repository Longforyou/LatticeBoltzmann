#! /usr/bin/env julia

# ======= Bounces

type Bounce{T <: Direction, V<:velocity_set._2D} <: Boundary
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}
end

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::Bounce{North})
    
    lbm.grid.f_prop[bound.rows, bound.cols, [4, 7, 8]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [2, 5, 6]]
end

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::Bounce{South})
    lbm.grid.f_prop[bound.rows, bound.cols, [2, 5, 6]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [4, 7, 8]]
end

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::Bounce{West})
    lbm.grid.f_prop[bound.rows, bound.cols, [3, 7, 6]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [1, 5, 8]]
end

function boundary(lbm::Lattice_Boltzmann_2D{Cells.D2Q9},
                  bound::Bounce{East})
    lbm.grid.f_prop[bound.rows, bound.cols, [1, 5, 8]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [3, 7, 6]]
end
