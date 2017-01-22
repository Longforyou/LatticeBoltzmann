#! /usr/bin/env Julia

using Plots, LatticeBoltzmann,
    LatticeBoltzmann._D2Q9
    D2Q9, North, South, West, East, Boundary, PeriodicPressure, Incompressible,
    LBM_Incompressible, LBM_Constants, Grid, Bounce, BGK, FullPeriodicStreaming,
    compute, get_pressure_pois, Collision

# Grid example
scale = 20
x = 5 * scale
y = 5 * scale
t =  4e4* x
write_inc = 5 #20 # After 50 Iter a file is created
U = 0.1 
H = 1.
L = 1. * H
nu_luft = 153.2e-7
mu_luft = 12.205e-6
rho_luft = 1.189
tau1 = sqrt(3/16) + 0.5

# Define all constants
consts =  LBM_Constants(U, H, nu_luft, rho_luft,
                       L, H, tau1)

# Create all objects for the LBM
grid = Grid_2D{D2Q9, Compressible}(consts, x, y, 9)


# Function to generate all objects for one simulation
function pois_compr(consts::LBM_Constants, grid::Grid_2D{D2Q9, Compressible})

    println("Constants: ", consts)
    bgk = BGK(consts)

    # Boundary conditions
    const top_bound = Bounce{North, D2Q9}(Array{Int64}(1:x), [1])
    const bottom_bound = Bounce{South, D2Q9}(Array{Int64}(1:x), [y])

    rho_out = rho_luft
    rho_in = rho_out + 3. * get_pressure_pois(consts)
    println("Rho_in: ", rho_in, "\n Rho_out: ", rho_out)
    peri_pres = PressurePeriodicStream_2D{West, D2Q9}(rho_in, rho_out,
                                             1, Array{Int64}(1:y),
                                                      x, Array{Int64}(1:y))
    full_stream = FullPeriodicStreaming_2D(grid)

    bounds = Array{Boundary, 1}([top_bound, bottom_bound])

    stream = Array{Streaming, 1}([peri_pres])
                             FullPeriodicStreaming,
                             BGK}(consts, grid, bounds)

    compute(grid, bgk, stream, bounds, "pois_", Array(1.:t), write_inc)

end

@time pois_compr(consts, grid)

# Analytical solution
y_vec = Array{Float64}((0:y-1))
uvec = get_velo_pois(Float64(x), Float64(y-1), consts, y_vec)

# Norm the results
y_vec ./= y
m_uvec = U / uvec[indmax(uvec)]
uvec ./= m_uvec

plot(y_vec, uvec, xlabel="y / y_max", ylabel="U_x/ U_max", label="analy",
     title="Analytical vs Numeric Pois")
plot!(y_vec,  grid.velocity[Int64(x), :, 2]./m_uvec, label="AuslassLBM Loesung")
plot!(y_vec,  grid.velocity[Int64(x/2), :, 2]./m_uvec, label="Mitte LBM Loesung")
plot!(y_vec,  grid.velocity[Int64(1), :, 2]./m_uvec, label="Einlass LBM Loesung")

savefig("pois_konvergenz.eps")

show()


