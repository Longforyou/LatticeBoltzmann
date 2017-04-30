#! usr/bin/julia

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
    ==


    
