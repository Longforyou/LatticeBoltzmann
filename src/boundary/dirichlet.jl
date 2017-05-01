#! /usr/bin/env julia

abstract Dirichlet <: Boundary

"""
    Pressure{T<:Direction, S<:SolType, V<:Velocity_Set} <: Boundary

Definition of a pressure/ density condition on the domain.
"""
immutable Pressure{T <: Direction, S <: SolType,
                   V <: Velocity_Set} <: Boundary
    rho::Float64
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}

end

"""
    compute_dirichlet{T<:Direction}(P, d2q9, f_i, post_ind, pre_ind, sign)

Kernel function for the computation of a dirichlet boundary condition.
"""
function compute_dirichlet{T<:Direction}(P::Pressure{T, NonEqBounce, D2Q9}, d2q9::D2Q9,
    f_i::Array{Float64, 1}, rho_ind::Array{Int64, 1},
    post_ind::Array{Int64, 1},
    pre_ind::Array{Int64, 1},
    sign::Bool)

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

"""
    boundary(grid, bound, d2q9)

Implementation of the interface function for the computation
of a dirichlet boundary condition. Dispatch 
"""
function boundary( grid::Grid_2D,
                   bound::Pressure{North, NonEqBounce,
                                   D2Q9}, d2q9::D2Q9)
    
    for row in bound.rows, col in bound.cols
        # Call the generitc function
        @fastmath @inbounds grid.f_prop[row, col, :] =
            compute_dirichlet(bound, grid.f_temp[row, col, :],
                  d2q9.dict[South], d2q9.dict[North],
                  [1, 2, 4], true )
    end 
    
end

function boundary( grid::Grid_2D, 
                   bound::Pressure{South, NonEqBounce,
                                   D2Q9}, d2q9::D2Q9)
    
    for row in bound.rows, col in bound.cols
        # Call the generic function
        @fastmath @inbounds grid.f_prop[row, col, :] =
            compute_dirichlet(bound, grid.f_temp[row, col, :],
                  d2q9.dict[North], d2q9.dict[South],
                  [1, 2, 4], true)
    end

end

function boundary( grid::Grid_2D,
                   bound::Pressure{West, NonEqBounce,
                                   D2Q9}, d2q9::D2Q9)


    for row in bound.rows, col in bound.cols
        # Call the generic function
        @fastmath @inbounds grid.f_prop[row, col, :] =
            compute_dirichlet(bound, grid.f_temp[row, col, :],
                  d2q9.dict[East], d2q9.dict[West],
                  [1, 3, 5], false)
    end

end

function boundary( grid::Grid_2D,
                   bound::Pressure{East, NonEqBounce,
                                   D2Q9}, d2q9::D2Q9)

    for row in bound.rows, col in bound.cols
        # Call the generic function
        @fastmath @inbounds grid.f_prop[row, col, :] =
            compute_dirichlet(bound, grid.f_temp[row, col, :],
                  d2q9.dict[West], d2q9.dict[East],
                  [1, 3, 5], true)
        
    end 
    
end
