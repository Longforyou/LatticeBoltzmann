#! /usr/bin/julia

# Define new abstract types for declaring 
abstract type absSolidComp end
abstract type absVirtualComp end
abstract type absBoundaryComp <: absVirtualComp end

include("lines.jl")
include("face.jl")
include("domain.jl")

export Line,
    Face,
    Domain,
    evalLineEq,
    linesParallel,
    linesIntersect,
    lineOnLine,
    pointOnLine,
    computeAngle,
    perpendicularLines
