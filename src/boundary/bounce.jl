#! /usr/bin/env julia

"""
    Bounce{T<:Direction, V<:Velocity} <: Boundary

Wrapper for the indices of a Bounceback condition on the domain.
"""
immutable Bounce{T <: Direction, V<:Velocity_Set} <: Boundary rows::Array{Int64, 1}
    cols::Array{Int64, 1}
end

function bounce_lattice!(lattices::Array{Lattice}, bounce::Bounce,
          to_arr::Array{Int64, 1}, from_arr::Array{Int64, 1})

    for row in bounce.rows, col in bounce.cols, i in 1:size(to_arr)[1]

        lattices[row, col].f_prop[to_arr[i]] =
            lattices[row, col].f_temp[from_arr[i]]::
    end
end

function boundary!(grid::Grid_2D,
                  bound::Bounce{North, D2Q9},
                  d2q9::D2Q9)
    
    # println("Pre\n", grid.lattices[5,1].f_prop)
    bounce_lattice!(grid.lattices, bound,
     d2q9.dict[South], d2q9.dict[North])
    
    # println("Post\n", grid.lattices[5,1].f_prop)
end

function boundary!(grid::Grid_2D,
                  bound::Bounce{South, D2Q9},
                  d2q9::D2Q9)
    bounce_lattice!(grid.lattices, bound,
     d2q9.dict[North], d2q9.dict[South])
end

function boundary!(grid::Grid_2D,
                  bound::Bounce{West, D2Q9},
                  d2q9::D2Q9)
                  
    grid.lattices[bound.rows, bound.cols].f_prop[dict[East]] =
      grid.lattices[bound.rows, bound.cols].f_temp[dict[West]]
end

function boundary!(grid::Grid_2D,
                  bound::Bounce{East, D2Q9},
                  d2q9::D2Q9)
    grid.lattices[bound.rows, bound.cols].f_prop[dict[West]] =
      grid.lattices[bound.rows, bound.cols].f_temp[dict[East]]
end
