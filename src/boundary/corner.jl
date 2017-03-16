#! /usr/bin/env julia

abstract Corner <: Boundary

"""
    Corner{V<:_2D}

Modelling of a corner for the 2D case. Assertion is done, since
there are only 4 possible cases. Pass a UInt8 out of [1, 2, 3, 4]
 corner 1 - 4
"""
immutable Corner_2D{V <: _2D} <: Corner

  row::Int64
  col::Int64
  quadrant::UInt8
  rho::Float64 # Specifies density/ pressure at the corner

<<<<<<< HEAD
=======
  Corner_2D(row, col, quad, rho) =
    (
     assert(quad in Array{UInt8}([1, 2, 3, 4])); # Check the quadrant
     new(row, col, quad, rho)
   )
>>>>>>> cd15aaf99060e3be94c3a55172a2b8dc94a8fa6c
end

# ============================================================
# ==== Corners
# ============================================================

"""
    boundary(grid, bound::Corner)

Implementation of the interface function for the computation of the
corner boundary condition. 
"""
function boundary(grid::Grid_2D, bound::Corner_2D{D2Q9})

  for row in bound.row, col in bound.col
    # Convenience pointer
    f_i = lbm.grid.f_temp[row, col, :]

    # Switch between the four cases
    if bound.quadrant == 0x01
        # f_3, f_7 & f_4 are known
        f_i[1] = f_i[3];  f_i[2] = f_i[4]
        f_i[5] = f_i[7];  f_i[6] = 0.; f_i[8] = 0.

    elseif bound.quadrant == 0x02
        # f_1, f_8 & f_4 are known
        f_i[3] = f_i[1];  f_i[2] = f_i[4]
        f_i[6] = f_i[8]; f_i[5] = 0.; f_i[7] = 0.

    elseif bound.quadrant == 0x03
        f_i[3] = f_i[1];  f_i[4] = f_i[2]
        f_i[7] = f_i[5];  f_i[5] = 0.; f_i[7] = 0.

    else
        f_i[1] = f_i[3];  f_i[4] = f_i[2]
        f_i[8] = f_i[6];  f_i[5] = 0.; f_i[7] = 0.

    end

    # Copy the changes over
    f_i[9] = bound.rho - sum(f_i[1:end-1])
    lbm.grid.f_temp[row, col, :] = f_i

   end
end
