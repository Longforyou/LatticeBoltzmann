#! /usr/bin/env julia

type LBM_Constants
 
  U::Float64   # Maximal Velocity
  Re::Float64  # Reynolds-Number
  char_L::Float64 # Characteristic length
  nu::Float64  # Viscocity (LBM)
  nu_F::Float64 # Viscosity (fluid)
  rho_F::Float64
  omega::Float64 # Relaxation factor
  phys_x::Float64 # Physical dimension of x
  phys_y::Float64
  
  LBM_Constants(_U, _char_L, nu_F, rho_F, _phyx, _phyy, _tau) =
  (
   _nu = (2. * _tau - 1.) / 6.;
   _Re = _char_L * _U / _nu;

   new(_U, _Re, _char_L, _nu, nu_F, rho_F,  1. / _tau,
      _phyx, _phyy)
  )
end

function get_rho(press::Float64)
  3 * press
end

function get_pressure_pois(consts::LBM_Constants)
  32. * consts.nu_F * consts.rho_F * consts.phys_x * consts.U / (consts.phys_y ^ 2) 
end

function get_velo_pois(L::Float64,
                       H::Float64,
                       consts::LBM_Constants, y::Array{Float64, 1})

    pres = get_pressure_pois(consts)
    pres ./ (2. * consts.nu * consts.rho_F *
              L) .* y .* (H - y)
              #-4. .* (consts.U / H^2) * .* y .* (H-y)

end
