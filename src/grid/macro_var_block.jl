#! usr/bin/julia

immutable SimpleMacroVariables2D{ F <: SingleDistFlow}

    density::SharedArray{ Float64, 2}
    velocity::SharedArray{ Float64, 3}

    SimpleMacroVariable2D(F::SingleDistFlow, length::Int64,
                          width::Int64) =
                              (
                                  density = SharedArray(Float64, (length, width));
                                  velocity = SharedArray(Float64, (length, width, 2));
                                  new(density, velocity)
                              )
end
