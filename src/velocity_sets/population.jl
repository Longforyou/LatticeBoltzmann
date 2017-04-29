#! usr/bin/julia



immutable LatticePopulation{F <: Flow, N <: Int64} <: Population_Set
    f_prop::SharedArray{Float64, N}
    f_eq::SharedArray{Float64, N}
    f_temp::Shared{Float64, N}
end
