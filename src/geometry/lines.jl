#! usr/env/julia

import Base: == 

immutable Line{T <: GeomProperty} <: Geometry
    startPoint::SVector
    koeff::SVector
    maxDistance::Float64

    Line(s::SVector, e::SVector, dist::Float64 = 0.25) =
        new(o, w - o, dist)
end

function evalLineEq(line::Line, t::Float64)
    return line.startPoint + t * line.koeff
end

function linesParallel(line1::Line, line2::Line)
    return line1.koeff == line2.koeff
end

function linesIntersect(line1::Line, line2::Line)
    return linesParallel(line1, line2) &&
        !(line1.startpoint != line2.startpoint)
end

function lineOnLine(line1::Line, line2::Line)
    if linesIntersect(line1, line2)
        b = -line1.startPoint -line2.startPoint
        A = hcat(line1.koeff, line2.koeff)
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
        return acos((line1.koeff * line2.koeff) /
                    (abs(line1.koeff) * abs(line2.koeff)))
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
