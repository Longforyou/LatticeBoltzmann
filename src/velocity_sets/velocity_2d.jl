#! /usr/bin/env julia

module _D2Q9

  using velocity_set
  abstract D2Q9 <: Particle

    # Some constant fields
    const c_x = Array{Float64,1}([0., 1., 0., -1., 0., 1., -1., -1., 1.])
    const c_y = Array{Float64,1}([0., 0., 1., 0., -1., 1., 1., -1., -1.])
    const w = Array{Float64,1}([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36])

export c_x, c_y, w, D2Q9

end # module _D2Q9
