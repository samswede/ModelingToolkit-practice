#using Pkg
#Pkg.add("ModelingToolkit")
#Pkg.add("OrdinaryDiffEq")
#Pkg.add("Plots")
#Pkg.add("PlotlyJS")
#Pkg.add("PlotlyBase")
#Pkg.add("PlotlyKaleido")

using ModelingToolkit
using LinearAlgebra
using OrdinaryDiffEq
using Plots
plotly()

function pendulumDouble!(du, u, p, t)
    x_1, dx_1, y_1, dy_1, T_1, x_2, dx_2, y_2, dy_2, T_2 = u
    g, L_1, L_2, m_1, m_2 = p

    du[1] = dx_1
    du[2] = T_1*x_1 + T_2*(x_2 - x_1) # unsure about sign of 2nd term
    du[3] = dy_1
    du[4] = T_1 * y_1 - g*m_1 + T_2*(y_2 - y_1) # unsure about sign of 2nd term
    du[5] = x_1^2 + y_1^2 - L_1^2

    du[6] = dx_2
    du[7] = T_2 * (x_2-x_1)
    du[8] = dy_2
    du[9] = T_2 * (y_2-y_1) - g*m_2
    du[10] = (x_2 - x_1)^2 + (y_2 - y_1)^2 - L_2^2

    return nothing
end

pendulum_fun! = ODEFunction(pendulum!, mass_matrix = Diagonal([1, 1, 1, 1, 0, 1, 1, 1, 1, 0]))

# Initial conditions (must be consistent?) If its a DAE with high index, the initial conditions must be consistent
u0 = [1.0, 0, 0, 0, 0, 2.0, 0, 0, 0, 0]
p = [9.8, 1, 1, 1, 1]
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
x_1 = sol[1, :]
y_1 = sol[3, :]

x_2 = sol[6, :]
y_2 = sol[8, :]

t = sol.t

# Generate the 3D plot
plot(t, x_1, y_1, xlabel="x", ylabel="y", zlabel="Time", title="3D Pendulum Motion", lw=2)

# Save the plot
#savefig("pendulum_3d_plot.png")

