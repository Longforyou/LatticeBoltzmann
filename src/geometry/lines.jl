#! usr/env/julia

import Base: ==, display

"""
    Line{T <: GeomProperty}

Type for the abstraction of a geometrical line. It
contains a start point and an maximal distance.
The interpolation is done with a coefficient computed with a passed end
point.
"""
mutable struct Line{T <: GeomProperty} <: Geometry
    startPoint::SVector
    coefficient::SVector
    maxDistance::Float64

    Line{T}(o, w, dist) where {T <: GeomProperty} =
        new(o, w - o, dist)
end

"""
    display(line)

Prints the values associated with the line.
"""
function display(line::Line)
    println("StartPoint: ", line.startPoint)
    println("Coefficient: ", line.coefficient)
    println("maximal Distance: ", line.maxDistance)
end

"""
    evalLineEq(line, t)

Computes the interpolated value for the passed line parameter t.
"""
function evalLineEq(line::Line, t::Float64)
    return line.startPoint + t * line.coefficient
end

"""
    linesParallel(line1, line2)

Returns true if and only if the lines are parallel to eachother
"""
function linesParallel(line1::Line, line2::Line)
    return line1.coefficient == line2.coefficient
end

"""
    linesIntersect(line1, line2)

Returns true if the lines are intersecting eachother.
"""
function linesIntersect(line1::Line, line2::Line)
    return linesParallel(line1, line2) &&
        !(line1.startpoint != line2.startpoint)
end

function lineOnLine(line1::Line, line2::Line)
    if linesIntersect(line1, line2)
        b = -line1.startPoint -line2.startPoint
        A = hcat(line1.coefficient, line2.coefficient)
        try
            w = inv(A'*A) * A'*A
        catch
            return false, 0
        end

        return true, w

    else
        return false, 0
    end
end

function pointOnLine(line::Line, p::SVector)
    temp_koef = (p - line.startPoint) /
        (line.endPoint - line.startPoint)
end

function computeAngle(line1::Line, line2::Line)
    if linesIntersect(line1, line2)
        return acos((line1.coefficient * line2.coefficient) /
                    (abs(line1.coefficient) * abs(line2.coefficient)))
    else
        return nothing
    end
end

function perpendicularLines(line1::Line, line2::Line)
    return computeAngle(line1, line2) == 90.
end

==(line1::Line, line2::Line) = (
    return lineOnLine(line1, line2)
)
