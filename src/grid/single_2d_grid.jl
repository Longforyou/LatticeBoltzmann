#! usr/bin/julia

immutable Simple2DGrid{ F <: SingleDistFlow} <: Grid

    length::Int64
    width::Int64
    directions::Int64
    latticePop::SingleLatticePopulation{ F, 3}
    macroVar::SimpleMacroVariable2D{ F }

    Simple2DGrid{F}(F::SingleDistFlow, length, width, directions) = (
        lattices = SingleLatticePopulation(F, length, width, directions);
        macroVar = SimpleMacroVariable2D( F, length, width );
        new(length, width, directions, lattices, macroVar)
    )

end
