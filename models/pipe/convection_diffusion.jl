# This one will formulate the equations based on forces
#import Pkg; Pkg.add("MethodOfLines"); Pkg.add("DomainSets")
using ModelingToolkit
using OrdinaryDiffEq
#using DifferentialEquations
using Plots
using MethodOfLines, DomainSets
plotly()


v= 0.05
C_in = 1
Daxial = 1e-3
L = 1
tf= 10

@parameters t, z

@variables C(..) # Find out why its not C(t, z), and how I could add other variables like v(t)
Dt = Differential(t)
Dz = Differential(z)
Dzz = Differential(z)^2


# equations
eq = Dt(C(t,z)) ~ Daxial*Dzz(C(t,z)) - v*Dz(C(t,z)) # for t > 0, z in (0, L)

#    eqs = [
#        Dt(C(t,z)) ~ Daxial*Dzz(C(t,z)) - v*Dz(C(t,z)), # for t > 0, z in (0, L)
#        C_in(t) ~ cos(t) # for t > 0
#        ]

# boundary conditions
boundary_conditions = [
    C(t, 0) ~ C_in, # for t > 0
    Dt(C(t, L)) ~ 0 # for t > 0
]

initial_conditions = [
    C(0, z) ~ 0 # for z in (0, L)
]

ic_bc = [
    # Boundary conditions
        # Inlet
    Dz(C(t, 0)) ~ v * (C(t, 0) - C_in) / Daxial,
    #C(t, 0) ~ C_in, # for t > 0 
        # Outlet
    Dz(C(t, L)) ~ 0, # for t > 0

    # Initial conditions
    C(0, z) ~ 1e-2 # for z in (0, L)
]

# Space and time domains
domains = [t ∈ (0.0,tf),
           z ∈ (0.0,L)]

@named pde_system = PDESystem(
    eq, #eqs
    ic_bc, #ic and bcs
    domains, #domain
    [t, z], #ivs
    [C(t, z)] #dvs
    )


# Extract the equations and unknowns
println("Equations")
println()
#println(equations(pde_system))
println()
#println("Unknowns")
#println()
#println(unknowns(pde_system))
println()
println("Parameters")
#println(parameters(pde_system))
println()
#println("Observed")
#println()
#println(observed(pde_system))


# Method of lines discretization
# Need a small dx here for accuracy
dz = 0.01
discretization = MOLFiniteDifference([z => dz], t)

println()
println("Made it past discretization")

# Convert the PDE problem into an ODE problem
prob = discretize(pde_system, discretization)
println()
println("Inspecting the type of prob")
println(typeof(prob)) # output: ODEProblem ...

# Note that it is at this point that I wish to add another ODE problem together with this one.

println()
println("Converted to ODE problem")

# Solve ODE problem
sol = solve(prob, Tsit5(), saveat=0.2)

println()
println("Solved ODE problem")

# Plot results and compare with exact solution
# Extract discrete values for z and t
discrete_z = sol[z]
discrete_t = sol[t]


# Plot 3D surface
surface(discrete_z, discrete_t, sol[C(t, z)], xlabel="z", ylabel="t", zlabel="C(t,z)", title="Concentration Profile")