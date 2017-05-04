#! usr/bin/julia

immutable D2Q9{F<:Flow} <: _2D
    
    w::Array{Float64, 1}
    dict::Dict{DataType, Array{Int64, 1}}

    D2Q9() = (
      # Some constant fields
      new(Array{Float64, 1}([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]),
          Dict{DataType, Array{Int64, 1}}([(North, [3, 6, 7]),
           (South, [5, 8, 9]), (West, [4, 7, 8]), (East, [2, 6, 9])]))
    )

 end
