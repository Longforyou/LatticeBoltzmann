#! usr/bin/julia

using RegionTrees

# A domain is a fixed data:
immutable Domain <: AbstractDomain
    blocks::Cell

end



