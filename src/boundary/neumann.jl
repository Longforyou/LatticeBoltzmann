#! /usr/bin/env julia


# ====== Neumann conditions
# The generic type T, allows for different equilibirum implementations
type Neumann{T<:Direction, S<:SolType, V <: Velocity_Set} <: Boundary
  vel::Float64
  rows::Array{Int64, 1}
  cols::Array{Int64, 1}

end

function compute_neumann(V::Velocity_Set._D2Q9.D2Q9,
                         f_i::Array{Float64, 1},
                         vel::Float64,
                         pre_ind::Array{Int64, 1},
                         post_ind::Array{Int64, 1},
                         vel_ind1::Array{Int64, 1}, sign::Bool)

  # Compute the new values
  if sign
    rho_i = (sum(f_i[vel_ind1]) + 2. *
             sum(f_i[pre_ind])) / (1. - vel)

    ru = rho_i * vel

    f_i[post_ind[1]] = f_i[pre_ind[1]] - (2./3.) * ru 
    f_i[post_ind[2]] = f_i[pre_ind[2]] - (1./6.) * ru +
      0.5 * (f_i[vel_ind1[2]] - f_i[vel_ind1[3]])
    f_i[post_ind[3]] = f_i[pre_ind[3]] - (1./6.) * ru +
      0.5 * (f_i[vel_ind1[3]] - f_i[vel_ind1[2]])
  
  else
    rho_i = (sum(f_i[vel_ind1]) + 2. *
             sum(f_i[pre_ind])) / (1. + vel)

    ru = rho_i * vel
    f_i[post_ind[1]] = f_i[pre_ind[1]] +
      (2./3.) * ru 
    f_i[post_ind[2]] = f_i[pre_ind[2]] + (1./6.) * ru -
      0.5 * (f_i[vel_ind1[2]] - f_i[vel_ind1[3]])
    f_i[post_ind[3]] = f_i[pre_ind[3]] + (1./6.) * ru -
      0.5 * (f_i[vel_ind1[3]] - f_i[vel_ind1[2]])

    end

  return f_i
end
