#! /usr/bin/env julia

module Cells

# definitions of types for referencing constant field etc.
abstract Particle
abstract Particle_2D <: Particle
abstract D2Q9 <:Particle_2D
  
module _2D

  using Cells
  abstract Particle_2D <:Cells.Particle

    module _D2Q9
    
    using Cells
    abstract D2Q9 <: Cells._2D.Particle_2D

    # Some constant fields
    const c_x = Array{Float64,1}([1., 0., -1., 0., 1., -1., -1., 1., 0.])
    const c_y = Array{Float64,1}([0., 1., 0., -1., 1., 1., -1., -1., 0.])
    const w = Array{Float64,1}([1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36,
                                       1/36, 4/9])

    end # module D2Q9

  end # module _2D

end # module Cells
