#using Pkg
#Pkg.add("ModelingToolkit")
#Pkg.add("OrdinaryDiffEq")
#Pkg.add("Plots")
#Pkg.add("PlotlyJS")

using ModelingToolkit
using LinearAlgebra
using OrdinaryDiffEq
using Plots
plotly()

function pendulum!(du, u, p, t)
    x, dx, y, dy, T = u
    g, L = p
    du[1] = dx
    du[2] = T * x
    du[3] = dy
    du[4] = T * y - g
    du[5] = x^2 + y^2 - L^2
    return nothing
end

pendulum_fun! = ODEFunction(pendulum!, mass_matrix = Diagonal([1, 1, 1, 1, 0]))

u0 = [1.0, 0, 0, 0, 0]
p = [9.8, 1]
tspan = (0, 10.0)

# Create specific instance of the problem, with initial conditions, time span and parameters
pendulum_prob = ODEProblem(pendulum_fun!, u0, tspan, p)

# Trace the problem to get the symbolic representation
traced_sys = modelingtoolkitize(pendulum_prob)

# Simplify the high-index system
# Use pantelides algorithm etc. to find an equivalent index-1 system
pendulum_sys = structural_simplify(dae_index_lowering(traced_sys))

# Reconstruct the problem from the simplified symbolic representation
prob = ODEProblem(pendulum_sys, [], tspan)

# Solve the problem
sol = solve(prob, Rodas5P(), abstol = 1e-8, reltol = 1e-8)

# Plot the solution
#plot(sol, idxs = unknowns(traced_sys))

# Extract the solution
x = sol[1, :]
y = sol[3, :]
t = sol.t

# Generate the 3D plot
plot(t, x, y, xlabel="x", ylabel="y", zlabel="Time", title="3D Pendulum Motion", lw=2)

# Save the plot
#savefig("pendulum_3d_plot.png")

