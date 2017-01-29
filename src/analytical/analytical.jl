#! /usr/bin/env Julia

"""
This file contains some analytical solution for a flow.
"""

function get_pressure_pois_1(consts::LBM_Constants)
  32. * consts.nu * 1. * consts.phys_x * consts.U / (consts.phys_y ^ 2) 
end

function get_pressure_pois_2(L::Float64, H::Float64, consts::LBM_Constants)
    3. * (L - 1.) * (8. * consts.nu * consts.U / H^2)
end

function get_velo_pois_1(L::Float64, H::Float64,
                       consts::LBM_Constants, y::Array{Float64, 1})

    pres = get_pressure_pois(consts)
    pres ./ (2. * consts.nu * 1. *  L) .* y .* (H - y)

end

function get_velo_pois_2(L::Float64, H::Float64,
                       consts::LBM_Constants, y::Array{Float64, 1})

       return -4. .* (consts.U / H^2) .* y .* (y - H)

end

export

    get_pressure_pois_1,
    get_pressure_pois_2,
    get_velo_pois_1,
    get_velo_pois_2
