#! /usr/bin/env Julia

using Plots, LatticeBoltzmann,
    LatticeBoltzmann._D2Q9

# Grid example
scale = 10
x = 20 * scale
y = 5 * scale
t =  80. * x
write_inc = 5 #:w20 # After 50 Iter a file is created
U = 0.1 
H = 4.
L = 1. * H
nu_luft = 153.2e-7
mu_luft = 12.205e-6
rho_luft = 1.189
tau1 = 3. * nu_luft + 0.5

# Create the velocity sets
d2q9 = D2Q9{Compressible}()

# Define all constants
consts =  LBM_Constants(U, H, nu_luft, rho_luft,
                       L, H, tau1)

# Create all objects for the LBM
grid = Grid_2D(consts, x, y, 9)


# Function to generate all objects for one simulation
function pois_compr(consts::LBM_Constants, grid::Grid_2D, d2q9::D2Q9{Compressible}, t::Float64)

    println("Constants: ", consts)
    bgk = BGK(consts)

    # Boundary conditions
    const top_bound = Bounce{North, D2Q9}(Array{Int64}(1:x), [1])
    const bottom_bound = Bounce{South, D2Q9}(Array{Int64}(1:x), [y])

    rho_out = 1.
    rho_in = rho_out + get_pressure_pois_2(Float64(x), Float64(y), consts)
    println("Rho_in: ", rho_in, "\n Rho_out: ", rho_out)

    peri_pres = PressurePeriodicStream_2D{West, D2Q9}(rho_in, rho_out,
                                             1, Array{Int64}(1:y),
                                                      x, Array{Int64}(1:y))

    full_stream = FullPeriodicStreaming_2D(grid)

    bounds = Array{Boundary, 1}([top_bound, bottom_bound])

    stream = Array{Streaming, 1}([peri_pres, full_stream])

    compute!(grid, d2q9, bgk, stream, bounds, "pois_air", t, write_inc)

end

@time pois_compr(consts, grid, d2q9, t)

# Analytical solution
y_vec = Array{Float64}(1:y) - 0.5
uvec = get_velo_pois_2(Float64(x), Float64(y), consts, y_vec)

# Norm the results
y_vec ./= y
m_uvec = uvec[indmax(uvec)]
uvec ./= m_uvec

plot(y_vec, uvec, xlabel="y / y_max", ylabel="U_x/ U_max", label="analytical",
     title="Analytische Lösung vs. LBM Lösung")
plot!(y_vec,  grid.velocity[Int64(x), :, 2]./m_uvec ./grid.density[Int64(x), :], label="AuslassLBM Loesung")
#plot!(y_vec,  grid.velocity[Int64(x/2), :, 2]./m_uvec, label="Mitte LBM Loesung")
plot!(y_vec,  grid.velocity[Int64(1), :, 2]./m_uvec ./grid.density[Int64(x), :], label="Einlass LBM Loesung")

savefig("pois_air.eps")

show()
