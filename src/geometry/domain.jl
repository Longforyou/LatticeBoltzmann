#! usr/bin/julia

# A domain is a fixed data:
struct Domain <: AbstractDomain
    blocks::Array{Float64, 2}
end



