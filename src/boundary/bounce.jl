#! /usr/bin/env julia

# ======= Bounces
type Bounce{T <: Direction, V<:velocity_set._D2Q9} <: Boundary
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}
end

function boundary(lbm::LBM{V <: velocity_set._D2Q9.D2Q9, F <: Flow,
                           S <: Streaming, C <: Collision},
                  bound::Bounce{North, V})
    
    lbm.grid.f_prop[bound.rows, bound.cols, [4, 7, 8]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [2, 5, 6]]
end

function boundary(lbm::LBM{V <: velocity_set._D2Q9.D2Q9, F <: Flow,
                           S <: Streaming, C <: Collision},
                  bound::Bounce{South})
    lbm.grid.f_prop[bound.rows, bound.cols, [2, 5, 6]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [4, 7, 8]]
end

function boundary(lbm::LBM{V <: velocity_set._D2Q9.D2Q9, F <: Flow,
                           S <: Streaming, C <: Collision},
                  bound::Bounce{West})
    lbm.grid.f_prop[bound.rows, bound.cols, [3, 7, 6]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [1, 5, 8]]
end

function boundary(lbm::LBM{V <: velocity_set._D2Q9.D2Q9, F <: Flow,
                           S <: Streaming, C <: Collision},
                  bound::Bounce{East})
    lbm.grid.f_prop[bound.rows, bound.cols, [1, 5, 8]] =
      lbm.grid.f_prop[bound.rows, bound.cols, [3, 7, 6]]
end
