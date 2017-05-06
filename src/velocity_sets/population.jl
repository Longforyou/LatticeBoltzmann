#! usr/bin/julia



immutable SingleLatticePopulation{ F <: SingleDistFlow,
                            N <: Int64 } <: SinglePopulationSet

    f_prop::SharedArray{ Float64, N }
    f_eq::SharedArray{ Float64, N }
    f_temp::SharedArray{ Float64, N }

    SingleLatticePopulation{ F, pops }(F::SingleDistFlow, length::Int64,
                                       width::Int64, pops::Int64) = (
                                           prop = SharedArray( Float64, (length, density, pops) );
                                           eq = SharedArray( Float64, (length, density, pops) );
                                           temp = SharedArray( Float64, (length, density, pops) );
                                           new(prop, eq, temp)
                                       )
end
