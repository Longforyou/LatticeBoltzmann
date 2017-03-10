#! /usr/bin/env julia

# ======= Bounces
type Bounce{T <: Direction, V<:Velocity_Set} <: Boundary
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}
end

function bounce_lattice!(lattices::Array{Lattice}, bounce::Bounce,
          to_arr::Array{Int64, 1}, from_arr::Array{Int64, 1})

    println("Bouncing")
    println("Pre\n", lattices[5].f_prop)
    for row in bounce.rows, col in bounce.cols
        @inbounds lattices[row, col].f_prop[to_arr] = lattices[row, col].f_temp[from_arr]
    end
    println("Post\n", lattices[5].f_prop)
end

function boundary!(grid::Grid_2D,
                  bound::Bounce{North, D2Q9},
                  d2q9::D2Q9)
    
    # grid.lattices[bound.rows, bound.cols].fprop[d2q9dict[South]] =
    #   grid.lattices[bound.rows, bound.cols].f_temp[d2q9dict[North]]
    bounce_lattice!(grid.lattices, bound,
     d2q9.dict[South], d2q9.dict[North])
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
