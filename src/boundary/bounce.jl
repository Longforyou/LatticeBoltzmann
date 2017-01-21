#! /usr/bin/env julia

# ======= Bounces
type Bounce{T <: Direction, V<:Velocity_Set} <: Boundary
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}
end

function boundary(lbm::LBM{ D2Q9, Flow, Streaming, Collision},
                  bound::Bounce{North, D2Q9})
    
    lbm.grid.f_prop[bound.rows, bound.cols, [4, 7, 8]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [2, 5, 6]]
end

function boundary(lbm::LBM{ D2Q9, Flow, Streaming, Collision},
                  bound::Bounce{South, D2Q9})
    lbm.grid.f_prop[bound.rows, bound.cols, [2, 5, 6]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [4, 7, 8]]
end

function boundary(lbm::LBM{ D2Q9, Flow, Streaming, Collision},
                  bound::Bounce{West, D2Q9})
    lbm.grid.f_prop[bound.rows, bound.cols, [3, 7, 6]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [1, 5, 8]]
end

function boundary(lbm::LBM{ D2Q9, Flow, Streaming, Collision},
                  bound::Bounce{East, D2Q9})
    lbm.grid.f_prop[bound.rows, bound.cols, [1, 5, 8]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [3, 7, 6]]
end
