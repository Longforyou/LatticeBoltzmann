#! /usr/bin/env Julia

using Plots, LatticeBoltzmann,
    LatticeBoltzmann._D2Q9

# Grid example
scale = 10
x = 20 * scale
y = 5 * scale
t =  40. * x
write_inc = 5 # After 50 Iter a file is created
U = Array(linspace(0.1, 0.3, 6))
H = 4.
L = 1. * H
nu_luft = 153.2e-7
mu_luft = 12.205e-6
rho_luft = 1.189

tau1 = sqrt(3/16) + 0.5

t_vec = 1.:t

# Create the velocity sets
d2q9 = D2Q9{Compressible}()

# Array for storing all different constants
const_vec = Array{LBM_Constants}(6)

# Define all constants
for i = 1:6

    const_vec[i] = LBM_Constants(U[i], H, nu_luft, rho_luft,
                        L, H, tau1)
end

println("Constants\n", const_vec)
# Create all objects for the LBM ( Grid needs only one arbitray LBM_Constants object for initialising the x, y values)
grid = Grid_2D(const_vec[1], x, y, 9)


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

    compute!(grid, d2q9, bgk, stream, bounds, "pois_cvelo", t)
    print(grid.velocity[end, :, 2] ./ grid.density[end, :])
    return grid.velocity[end, :, 2] ./ grid.density[end, :]

end



# Analytisch solution
y_vec = Array{Float64}(1:y) - 0.5

# Plot all 6 values
# ===== 1
uvec = get_velo_pois_2(Float64(x), Float64(y), const_vec[1], y_vec)
# Norm the results
y_vec ./= y
m_uvec = const_vec[1].U
uvec ./= m_uvec
p1 = plot(xlabel="y", yaxis="U_x", title=string("LBM_Lösung", const_vec[1].U))
plot!(p1, y_vec, uvec, label="analytisch")
plot!(p1, y_vec, pois_compr(const_vec[1], grid, d2q9, t) ./ m_uvec, label=string("U: ", const_vec[1].U))


y_vec = Array{Float64}(1:y) - 0.5
uvec = get_velo_pois_2(Float64(x), Float64(y), const_vec[2], y_vec)
# Norm the results
y_vec ./= y
m_uvec = const_vec[2].U
uvec ./= m_uvec
p2 = plot(xlabel="y", yaxis="U_x", title=string("U=", const_vec[2].U))
plot!(p2, y_vec, uvec, label="analytisch")
plot!(p2, y_vec, pois_compr(const_vec[2], grid, d2q9, t) ./ m_uvec, label=string("U: ", const_vec[2].U))

y_vec = Array{Float64}(1:y) - 0.5
uvec = get_velo_pois_2(Float64(x), Float64(y), const_vec[3], y_vec)
# Norm the results
y_vec ./= y
m_uvec = const_vec[3].U
uvec ./= m_uvec
p3 = plot(xlabel="y", yaxis="U_x", title=string("U=", const_vec[3].U))
plot!(p3, y_vec, uvec, label="analytisch")
plot!(p3, y_vec, pois_compr(const_vec[3], grid, d2q9, t) ./ m_uvec, label=string("U: ", const_vec[3].U))

y_vec = Array{Float64}(1:y) - 0.5
uvec = get_velo_pois_2(Float64(x), Float64(y), const_vec[4], y_vec)
# Norm the results
y_vec ./= y
m_uvec = const_vec[4].U
uvec ./= m_uvec
p4 = plot(xlabel="y", yaxis="U_x", title=string("U=", const_vec[4].U))
plot!(p4, y_vec, uvec, label="analytisch")
plot!(p4, y_vec, pois_compr(const_vec[4], grid, d2q9, t) ./ m_uvec, label=string("U: ", const_vec[4].U))

y_vec = Array{Float64}(1:y) - 0.5
uvec = get_velo_pois_2(Float64(x), Float64(y), const_vec[5], y_vec)
# Norm the results
y_vec ./= y
m_uvec = const_vec[5].U
uvec ./= m_uvec
p5 = plot(xlabel="y", yaxis="U_x", title=string("U=", const_vec[5].U))
plot!(p5, y_vec, uvec, label="analytisch")
plot!(p5, y_vec, pois_compr(const_vec[5], grid, d2q9, t) ./ m_uvec, label=string("U: ", const_vec[5].U))

y_vec = Array{Float64}(1:y) - 0.5
uvec = get_velo_pois_2(Float64(x), Float64(y), const_vec[6], y_vec)
# Norm the results
y_vec ./= y
m_uvec = const_vec[6].U
uvec ./= m_uvec
p6 = plot(xlabel="y", yaxis="U_x", title=string("U=", const_vec[6].U))
plot!(p6, y_vec, uvec, label="analytisch")
plot!(p6, y_vec, pois_compr(const_vec[6], grid, d2q9, t) ./ m_uvec, label=string("U: ", const_vec[6].U))

plot(p1, p2, p3, p4, p5, p6)
# Saving..
savefig("pois_velo_diff.eps")

show()
