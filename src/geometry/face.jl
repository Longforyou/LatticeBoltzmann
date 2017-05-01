#! usr/bin/julia

immutable Face{T <: GeomProperty} <: Geometry
    lines::Array{Line}
    angles::Array{Float64}

    Face(lines::Array{Line}) = (
        linelength = length(lines);
        assert(linelength >= 4);
        angles = zeros(linelength);
        for i = 2:length(lines);
           angles[i-1] = computeAngle(lines[i-1], line[i]);
        end;
        angles[end] = computeAngle(lines[end], lines[1]);
        new(lines, angles)
    )
end


