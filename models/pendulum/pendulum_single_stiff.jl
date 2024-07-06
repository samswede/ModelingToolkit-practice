# This one will formulate the equations based on forces

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using LinearAlgebra
using OrdinaryDiffEq
using DifferentialEquations
using Plots
plotly()

@mtkmodel Pendulum begin
    @parameters begin
        xa=0
        ya=0
        m=1.0
        L0=1.0
        k=1e6
        g=9.8
        F_g = m*g
    end
    @variables begin
        x(t)
        y(t)

        v_x(t)
        v_y(t)

        L(t)

        T(t)
        T_x(t)
        T_y(t)
    end
    @equations begin
        # connectors
        #port.xa ~ xa
        #port.ya ~ ya


        # Component Forces
        T_x ~ T*((x - xa)/L)
        T_y ~ T*((y - ya)/L)
        #F_g ~ m*g

        D(x) ~ v_x
        D(y) ~ v_y
        # physics
        m*D(v_x)~ -(T_x)
        m*D(v_y) ~ -T_y - F_g

        (x-xa)^2 + (y-ya)^2 ~ L^2 # Constraint equation

        T ~ k*(L - L0) # Spring force (compliance method)

    end
end

#when using @mtkbuild then structural_simplify is automatically called and we therefore cannot see the unsimplify system. Replace @mtkbuild with @named to generate an ODESystem without applying structural_simplify.
@mtkbuild pendulum = Pendulum() # automatically applies structural_simplify
     # automatically applies structural_simplify

# Extract the equations and unknowns
println("Equations")
println()
println(equations(pendulum))
println()
println("Unknowns")
println()
println(unknowns(pendulum))
println()
println("Parameters")
println(parameters(pendulum))
println()
println("Observed")
println()
println(observed(pendulum))


function set_initial_conditions!(pendulum, θ, L, xa, ya)
    # Calculate initial positions based on the angle θ and length L
    x_init = xa + L * sin(θ)
    y_init = ya - L * cos(θ)
    
    # Initial conditions
    initial_conditions = [
        pendulum.x => x_init,
        pendulum.y => y_init,
        pendulum.L => L,
        pendulum.v_x => 0.0,
        pendulum.v_y => 0.0
    ]
    
    # Initial guesses for the remaining variables
    initial_guesses = [
        pendulum.T => 0, 
        pendulum.T_x => 0, 
        pendulum.T_y => 0.0
    ]
    
    return initial_conditions, initial_guesses
end

# Example parameters
θ = pi/4  # 45 degrees
length = 10.0
location_x = 100.0
location_y = 100.0

# Set initial conditions using the function
initial_conditions, initial_guesses = set_initial_conditions!(pendulum, θ, length, location_x, location_y)


tspan = (0, 10)

pendulum.xa = location_x
pendulum.ya = location_y
pendulum.L0 = length

# Reconstruct the problem from the simplified symbolic representation
prob = ODEProblem(pendulum, initial_conditions, tspan, guesses = initial_guesses)
# Note that you need to use pendulum.x etc instead of x, y etc

dt=1e-6
# Solve the problem
sol = solve(prob, Rodas5P(), initializealg=ShampineCollocationInit(dt))

# Plot the solution
#plot(sol, idxs = unknowns(traced_sys))

# Extract the solution
x = sol[pendulum.x, :]
y = sol[pendulum.y, :]


# Generate the 3D plot
plot(sol.t, x, y, xlabel="Time", ylabel="y", zlabel="x", title="3D Pendulum Motion", lw=2)


