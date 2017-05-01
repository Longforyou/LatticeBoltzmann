#! usr/bin/julia

# Functions to test
function lineEvalTest()
    line = Line{Solid}(SVector(1., 0., 0.),
                       SVector(3., 0., 0.))

    t::Float64 = 4.
    testVal = evalLineEq(line, t)

    return testVal == SVector(8., 0., 0.)
end

function lineParallelTest()
    line = Line{Solid}(SVector(1., 0., 0.),
                       SVector(3., 0., 0.))
    line2 = Line{Solid}(SVector(2., 0., 0.),
                       SVector(4., 0., 0.))

    return linesParallel(line, line2)

end

function lineIntersectTest()
    line = Line{Solid}(SVector(1., 0., 0.),
                       SVector(3., 0., 0.))
    line2 = Line{Solid}(SVector(4., 0., 0.),
                       SVector(4., 1., 0.))

    return linesParallel(line, line2)
end

function lineOnLineTest()
    line = Line{Solid}(SVector(1., 0., 0.),
                       SVector(3., 0., 0.))
    line2 = Line{Solid}(SVector(2., 0., 0.),
                       SVector(4., 0., 0.))

    check, p = lineOnLine(line, line2) 

    return check

end

function computeAngleTest()
    line = Line{Solid}(SVector(1., 0., 0.),
                       SVector(3., 0., 0.))
    line2 = Line{Solid}(SVector(2., 0., 0.),
                       SVector(4., 0., 0.))

    check = computeAngle(line, line2) == 0.
end

function perpendicularLineTest()
    line = Line{Solid}(SVector(1., 0., 0.),
                       SVector(3., 0., 0.))
    line2 = Line{Solid}(SVector(2., 0., 0.),
                        SVector(4., 1., 0.))

    return perpendicularLines(line, line2)
end

@testset "Geometry Test" begin
    @testset "Line Eval Test" begin
        @test lineEvalTest()
    end

    @testset "Line Parallel Test" begin
        @test lineParallelTest()
    end

    @testset "Intersect Test" begin
        @test intersectTest()
    end

    @testset "Line On Line Test" begin
        @test lineOnLineTest()
    end

    @testset "Angle Test" begin
        @test computeAngleTest()
    end

    @testset "Perpendicular Test" begin
        @test perpendicularLineTest()
    end
end
