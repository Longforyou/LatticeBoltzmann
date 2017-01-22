#! /usr/bin/env julia

"""
This file contains some analytical solution for a flow.
"""

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

export
    get_pressure_pois,
    get_velo_pois
