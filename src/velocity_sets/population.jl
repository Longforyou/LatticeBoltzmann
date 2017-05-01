#! usr/bin/julia



immutable LatticePopulation{F <: Flow} <: Population_Set
    f_prop::SharedArray{Lattice, N}
    f_eq::SharedArray{Lattice, N}
    f_temp::SharedArray{Lattice, N}
end
