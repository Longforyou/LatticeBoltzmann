#! /usr/bin/env julia


# ======= Pressure Dirichlet Conditions
abstract Dirichlet <: Boundary

type Pressure{T <: Direction, S <: SolType, V<:Velocity_Set} <: Boundary
  rho::Float64
  rows::Array{Int64, 1}
  cols::Array{Int64, 1}
   

end

function compute_dirichlet(V::D2Q9, f_i::Array{Float64, 1}, rho::Float64,
                           pre_ind::Array{Int64, 1}, 
                           post_ind::Array{Int64, 1},
                           rho_ind::Array{Int64, 1}, sign::Bool)

    ru = rho - abs(sum(f_i[rho_ind]) + 2. *
                   sum(f_i[pre_ind]))

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
