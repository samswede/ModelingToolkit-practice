using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using LinearAlgebra
using OrdinaryDiffEq
using DifferentialEquations
using Plots
plotly()


@mtkmodel ImmovablePin begin
    @parameters begin
        a=0.0
        b=0.0
    end
    @variables begin
        T_x(t)
        T_y(t)
    end
    @components begin
        port = PendulumPort()
    end
    @equations begin

        #v_y ~ 0
        #v_x ~ 0

        # Connector
        port.T_x ~ -T_x # The force is opposite to the tension in the pendulum
        port.T_y ~ -T_y # The force is opposite to the tension in the pendulum
        port.v_xa ~ 0
        port.v_ya ~ 0
        
    end
end


@connector PendulumPort begin
    # Across variables
    v_xa(t) 
    v_ya(t) 

    # Through variables (flow variables)
    T_x(t), [connect = Flow] # Tension in the pendulum
    T_y(t), [connect = Flow] # Tension in the pendulum

end

@connector MechanicalPort begin
    # Across variables
    v_x(t) # Velocity in the x direction
    v_y(t) # Velocity in the y direction

    # Through variables (flow variables)
    Fx(t), [connect = Flow] # Force in the x direction
    Fy(t), [connect = Flow] # Force in the y direction

end


@mtkmodel Pendulum begin
    @parameters begin
        m=1.0
        g=9.8

        L0=1.0
        k=1e6

        F_g = m*g
    end
    @variables begin
        xa(t)
        ya(t)
        v_xa(t)
        v_ya(t)

        x(t)
        y(t)
        v_x(t)
        v_y(t)

        L(t)

        T(t)
        T_x(t)
        T_y(t)
    end
    @components begin
        port_a = PendulumPort()
        #port_b = MechanicalPort()
    end
    @equations begin

        D(x) ~ v_x
        D(y) ~ v_y

        D(xa) ~ v_xa
        D(ya) ~ v_ya

        # Component Forces
        T_x ~ T*((x - xa)/L)
        T_y ~ T*((y - ya)/L)
        #F_g ~ m*g


        # physics
        m*D(v_x)~ -(T_x)
        m*D(v_y) ~ -T_y - F_g

        (x-xa)^2 + (y-ya)^2 ~ L^2 # Constraint equation

        T ~ k*(L - L0) # Spring force (compliance method)


        # connectors
        port_a.v_xa ~ v_xa
        port_a.v_ya ~ v_ya
        port_a.T_x ~ -T_x # The force is opposite to the tension in the pendulum
        port_a.T_y ~ -T_y # The force is opposite to the tension in the pendulum

    end
end

@mtkmodel System begin
    @parameters begin
        a=0.0
        b=0.0
        
    end
    @components begin
        immovable_pin = ImmovablePin()
        pendulum = Pendulum()
    end
    @equations begin
        connect(immovable_pin.port, pendulum.port_a)

    end

end

@mtkbuild system = System() 

# Extract the equations and unknowns
println("Equations")
println()
println(equations(system))
println()
println("Unknowns")
println()
println(unknowns(system))
println()
println("Parameters")
println(parameters(system))
println()
println("Observed")
println()
println(observed(system))



function set_initial_conditions!(system, θ, L, xa, ya)
    # Calculate initial positions based on the angle θ and length L
    x_init = xa + L * sin(θ)
    y_init = ya - L * cos(θ)
    
    # Initial conditions
    initial_conditions = [
        system.pendulum.x => x_init,
        system.pendulum.y => y_init,
        
        system.pendulum.L => L,
        system.pendulum.v_x => 0.0,
        system.pendulum.v_y => 0.0,
        system.pendulum.xa => 0

    ]
    

    # Initial guesses for the remaining variables
    initial_guesses = [
        
        system.pendulum.T => 0, 
        system.pendulum.T_x => 0, 
        system.pendulum.T_y => 0.0,

        system.pendulum.ya => 0,
        system.pendulum.v_xa => 0.0,
        system.pendulum.v_ya => 0.0
        
    ]
    
    return initial_conditions, initial_guesses
end

# Example parameters
θ = pi/4  # 45 degrees
length = 10.0
x0 = 100.0
y0 = 100.0
tspan = (0, 10)

# Set initial conditions using the function
initial_conditions, initial_guesses = set_initial_conditions!(system, θ, length, x0, y0)


# Reconstruct the problem from the simplified symbolic representation
prob = ODEProblem(system, initial_conditions, tspan, guesses = initial_guesses)
# Note that you need to use pendulum.x etc instead of x, y etc

dt=1e-6
# Solve the problem
sol = solve(prob, Rodas5P(), initializealg=ShampineCollocationInit(dt))

# Plot the solution
#plot(sol, idxs = unknowns(traced_sys))

# Extract the solution
x = sol[system.pendulum.x, :]
y = sol[system.pendulum.y, :]


# Generate the 3D plot
plot(sol.t, x, y, xlabel="x", ylabel="y", zlabel="Time", title="3D system Motion", lw=2)

