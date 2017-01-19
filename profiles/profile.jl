#! /usr/bin/env julia

# NOTES
# Runs this script inside a julia REPl, with
# include("profile.jl")
# After the script has finished it's computation, write
# ProfileView.view()
# to see the Firediagramm of the function
#
using ProfileView

include("setup.jl")  # Load all modules

# Grid example
x = 100
y = 50
t = 1000
write_inc = 20 # After 50 Iter a file is created

function profile_func(n::Int)
  grid = Grid{Cells._2D._D2Q9.D2Q9}(x, y, 9)
  consts = LBM_Constants(0.1/y, x, 2.0, 1.0, sqrt(3/16) + 0.5)

  #top_boun = Boundary_Conditions.BounceCondition(Array(1:x), [1], "North")
  #bottom_boun = Boundary_Conditions.BounceCondition(Array(1:x), [y], "South")
  top_boun = Boundary_Conditions.BounceNorth(Array{Int64}(1:x), [1])
  bottom_boun = Boundary_Conditions.BounceSouth(Array{Int64}(1:x), [y])

  side_vent = Boundary_Conditions.PeriodicPressureCondition(1, Array{Int64}(1:y),
                                                            x, Array{Int64}(1:y), 1., 0)
  get_periodic_pressure(side_vent, consts)

  bounds = Array{Boundary_Conditions.Boundary, 1}([side_vent, bottom_boun, top_boun])

  lbm = Lattice_Boltzmann_2D{Cells._2D._D2Q9.D2Q9}(consts, grid, bounds)

  compute(lbm, "impl_test", Array(1.0:t), write_inc)
end

profile_func(1)
Profile.clear()
@profile profile_func(10)
#ProfileView.view()
#ProfileView.svgwrite("profile_results.svg")
