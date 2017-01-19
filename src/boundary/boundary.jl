#! /usr/bin/env julia

#= 
Basic idea is provide indices during the computations.
This provides other function to directly compute the data

=#

# Abstract type for the dispatch
abstract Boundary
abstract Dirichlet <: Boundary
abstract SolType
abstract Eq <: SolType
abstract NonEq <: SolType
abstract NonEqBounce <: SolType
abstract BounceCondition <: Boundary
abstract Direction
abstract North <: Direction
abstract South <: Direction
abstract West <: Direction
abstract East <: Direction
#abstract Pressure{SolType} <: Boundary

# ====== Neumann conditions
# The generic type T, allows for different equilibirum implementations
type Neumann{T<:Direction, S<:SolType} <: Boundary# <: Neumann{T<:SolType}
  vel::Float64
  rows::Array{Int64, 1}
  cols::Array{Int64, 1}

end

function compute_neumann(f_i::Array{Float64, 1}, vel::Float64,
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

# ======= Pressure Dirichlet Conditions

# The generic type T, allows for different equilibirum implementations

type Pressure{T <: Direction, S <: SolType} <: Boundary
  rho::Float64
  rows::Array{Int64, 1}
  cols::Array{Int64, 1}

end

function compute_dirichlet(f_i::Array{Float64, 1}, rho::Float64,
                           pre_ind::Array{Int64, 1}, 
                           post_ind::Array{Int64, 1},
                           rho_ind::Array{Int64, 1}, sign::Bool)

    ru = rho - abs(sum(f_i[rho_ind]) + 2. *
                   sum(f_i[pre_ind]))

    # println("ru: ", ru)
    # println("Pre Mod: ", f_i)

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
    # println("Post Mod: ", f_i)

    return f_i
end

type PeriodicPressure{T <: Direction} <: Boundary
    rho_inlet::Float64
    rho_outlet::Float64
    inlet_row::Int64
    inlet_col::Array{Int64, 1}
    outlet_row::Int64
    outlet_col::Array{Int64, 1}

end

# ======= Bounces

type Bounce{T <: Direction} <: Boundary
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}
end

# Open Boundary
type OpenBound{T <: Direction} <: Boundary
    rows::Array{Int64, 1}
    cols::Array{Int64, 1}
  end


# The corners are modelling different quadrant for the grid. By the type of the
# corner 1 - 4
type Corner <: Boundary

  row::Int64
  col::Int64
  quadrant::UInt8
  rho::Float64 # Specifies density/ pressure at the corner

  Corner(row, col, quad, rho) =
    (
     assert(quad in Array{UInt8}([1, 2, 3, 4])); # Check the quadrant
     new(row, col, quad, rho)
    )
end



