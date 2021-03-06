#! /usr/bin/env Julia

using Plots, LatticeBoltzmann,
    LatticeBoltzmann._D2Q9

# Grid example
scale = 10
x = 20 * scale
y = 5 * scale
t =  40. * x
write_inc = 5 # After 50 Iter a file is created
U = 0.1 
H = 4.
L = 1. * H
nu_luft = 153.2e-7
mu_luft = 12.205e-6
rho_luft = 1.189

# Define all relaxation terms
tau_vec = Array(linspace(0.1, 2, 5))
tau1 = insert!(tau_vec, 4, sqrt(3/16) + 0.5)

# Define the array for storing all time histories
t_vec = 1.:t

# Create the velocity sets
d2q9 = D2Q9{Compressible}()

# Array for storing all different constants
const_vec = Array{LBM_Constants}(6)

# Array for storing the errors

# Define all constants
for i = 1:6

    const_vec[i] =  LBM_Constants(U, H, nu_luft, rho_luft,
                        L, H, tau_vec[i])
end

# Create all objects for the LBM ( Grid needs only one arbitray LBM_Constants object for initialising the x, y values)
grid = Grid_2D(const_vec[1], x, y, 9)


# Function to generate all objects for one simulation
function pois_compr(consts::LBM_Constants, grid::Grid_2D, d2q9::D2Q9{Compressible}, uvec::Array{Float64, 1}, t::Float64)

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

    return compute!(grid, d2q9, bgk, stream, bounds, "pois_conv", t, uvec, 2)

end



# Analytical solution
y_vec = Array{Float64}(1:y) - 0.5

plot(xlabel="t", yaxis=(:log, "L^2 Error"), title="Fehler im Verlauf der Zeit")

for const_i in const_vec

    uvec = get_velo_pois_2(Float64(x), Float64(y), const_i, y_vec)
    plot!(t_vec, pois_compr(const_i, grid, d2q9, uvec, t), 
        label=string("Tau= ", const_i.tau, " Nu= ", const_i.nu))
end

savefig("pois_tau_konvergenz.eps")

show()
