#! /usr/bin/env julia


# ======= Pressure Dirichlet Conditions
abstract Dirichlet <: Boundary

immutable Pressure{T <: Direction, S <: SolType,
                   V <: Velocity_Set} <: Boundary
    rho::Float64
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}

end

function (P::Pressure{Direction, NonEqBounce, D2Q9})(
    f_i::Array{Float64, 1}, rho_ind::Array{Int64, 1}, pre_ind::Array{Int64, 1},
    post_ind::Array{Int64, 1}, sign::Bool)

    ru = -1. + (sum(f_i[rho_ind]) + 2. * sum(f_i[pre_ind])) / P.rho

    if sign
        f_i[post_ind[1]] = f_i[pre_ind[1]] - 
            (2./3.) * ru 
        f_i[post_ind[2]] = f_i[pre_ind[2]] - 
            (1./6.) * ru + 0.5 * (f_i[rho_ind[2]] - f_i[rho_ind[3]])
        f_i[post_ind[3]] = f_i[pre_ind[3]] - 
            (1./6.) * ru + 0.5 * (f_i[rho_ind[3]] - f_i[rho_ind[2]])
    else

        f_i[post_ind[1]] = f_i[pre_ind[1]] + 
            (2./3.) * ru 
        f_i[post_ind[2]] = f_i[pre_ind[2]] + 
            (1./6.) * ru + 0.5 * (f_i[rho_ind[2]] - f_i[rho_ind[3]])
        f_i[post_ind[3]] = f_i[pre_ind[3]] + 
            (1./6.) * ru + 0.5 * (f_i[rho_ind[3]] - f_i[rho_ind[2]])
    end

    return f_i
end

# ===========================================================
# === Dirichlet(Pressure) solution schemes
# ===========================================================
function boundary( grid::Grid_2D,
                   bound::Pressure{North, NonEqBounce,
                                   D2Q9})
    
    for row in bound.rows, col in bound.cols
        # Call the generitc function
        lbm.grid.f_temp[row, col, :] =
            bound(lbm.grid.f_temp[row, col, :],
                  [5, 8, 9], [3, 6, 7],
                  [1, 2, 4], true )
    end 
    
end

function boundary( grid::Grid_2D, 
                   bound::Pressure{South, NonEqBounce,
                                   D2Q9})
    
    for row in bound.rows, col in bound.cols
        # Call the generic function
        lbm.grid.f_temp[row, col, :] =
            bound(lbm.grid.f_temp[row, col, :],
                  [3, 6, 7], [5, 8, 9],
                  [1, 2, 4], true)
    end

end

function boundary( grid::Grid_2D,
                   bound::Pressure{West, NonEqBounce,
                                   D2Q9})


    for row in bound.rows, col in bound.cols
        # Call the generic function
        lbm.grid.f_temp[row, col, :] =
            bound(lbm.grid.f_temp[row, col, :],
                  [4, 8, 7], [2, 6, 9],
                  [1, 3, 5], false)
    end

end

function boundary( grid::Grid_2D,
                   bound::Pressure{East, NonEqBounce,
                                   D2Q9})

    for row in bound.rows, col in bound.cols
        # Call the generic function
        lbm.grid.f_temp[row, col, :] =
            bound(lbm.grid.f_temp[row, col, :],
                  [2, 6, 9], [4, 8, 7],
                  [1, 3, 5], true)
        
    end 
    
end
