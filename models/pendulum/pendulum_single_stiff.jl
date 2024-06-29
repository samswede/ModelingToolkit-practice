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
        x0=0
        y0=0
        m=1.0
        g=9.8

        L0=1.0
        k=1e6

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
        #port.x0 ~ x0
        #port.y0 ~ y0

        # Component Forces
        T_x ~ T*((x - x0)/L)
        T_y ~ T*((y - y0)/L)
        #F_g ~ m*g

        D(x) ~ v_x
        D(y) ~ v_y
        # physics
        m*D(v_x)~ -(T_x)
        m*D(v_y) ~ -T_y - F_g

        x^2 + y^2 ~ L^2 # Constraint equation

        T ~ k*(L - L0) # Spring force (compliance method)

    end
end

#when using @mtkbuild then structural_simplify is automatically called and we therefore cannot see the unsimplify system. Replace @mtkbuild with @named to generate an ODESystem without applying structural_simplify.
@mtkbuild pendulum = Pendulum() # automatically applies structural_simplify

#@named pendulum = ODESystem(Pendulum()) # does not apply structural_simplify
# perform index reduction and simplify the system
#pendulum = structural_simplify(dae_index_lowering(pendulum))

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


tspan = (0, 10)

# The following ICs must be satisfied exactly
initial_conditions = [
    pendulum.x => 1.0, 
    pendulum.y => 0.0, 
    pendulum.L => 1.0, 
    pendulum.v_x => 0.0, 
    pendulum.v_y => 0.0
    ]

# The remaining variables must be given initial guesses
initial_guesses = [
    pendulum.T => 0, 
    pendulum.T_x => 0, 
    pendulum.T_y => 0.0
    ]

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
plot(sol.t, x, y, xlabel="x", ylabel="y", zlabel="Time", title="3D Pendulum Motion", lw=2)


