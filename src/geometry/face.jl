#! usr/bin/julia

struct Face{T <: GeomProperty} <: Geometry
    lines::Array{Line}
    angles::Array{Float64}

    Face{T}(lines::Array{Line}) where {T <: GeomProperty}= (
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


