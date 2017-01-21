
# The corners are modelling different quadrant for the grid. By the type of the
# corner 1 - 4
type Corner{V <: _2D} <: Boundary

  row::Int64
  col::Int64
  quadrant::UInt8
  rho::Float64 # Specifies density/ pressure at the corner

  # Corner{V}(row, col, quad, rho) =
  #   (
  #    assert(quad in Array{UInt8}([1, 2, 3, 4])); # Check the quadrant
  #    new(row, col, quad, rho)
  #  )
end

# ============================================================
# ==== Corners
# ============================================================
# TODO Change the indices
function boundary(lbm::LBM{D2Q9, Flow, Streaming, Collision},
                  bound::Corner{D2Q9})

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
