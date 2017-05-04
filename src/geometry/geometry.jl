#! /usr/bin/julia

using GeometricTypes

# Define new abstract types for declaring 
abstract absSolidComp
abstract absVirtualComp
abstract absBoundaryComp <: absVirtualComp

include("lines.jl")
include("face.jl")
include("domain.jl")

export
    Line,
    Face,
    Domain,

    # Functions
    evalLineEq,
    linesParallel,
    linesIntersect,
    lineOnLine,
    pointOnLine,
    computeAngle,
    perpendicularLines,


    
