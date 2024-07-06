using MethodOfLines, ModelingToolkit, DomainSets, OrdinaryDiffEq, Plots

function test(L::Float64, T::Float64)
    @parameters τ, x
    @variables q(..), c(..)  

    ∂τ = Differential(τ) 
    ∂x = Differential(x)
    ∂xx = Differential(x)^2
    
    ϵ = 0.4436 
    b = 0.0258
    k_f = 0.89
    D = 2.69e-4
    u = 0.0424
    qm = 76.85
    c0 = 8.8

    eq = [
        0 ~ D * ∂xx(c(τ, x)) - ∂τ(c(τ, x)) + (1 - ϵ) / ϵ * (-∂τ(q(τ, x))) + u * ∂x(c(τ, x)),
        0 ~ -∂τ(q(τ, x)) + k_f * (q(τ, x) + qm * b * c(τ, x) / (1 + b * c(τ, x)))
    ]

    domain = [x ∈ Interval(0.0, L),
              τ ∈ Interval(0.0, T)
    ]

    ic_bc = [
        # initial conditions
        q(T, x) ~ 0.0,
        c(T, x) ~ 0.0,
        # boundary conditions
        ∂x(c(τ, 0)) ~ u * (c(τ, 0) - c0) / D,
        ∂x(c(τ, L)) ~ 0.0
    ]

    @named pde_system = PDESystem(
        eq, ic_bc, domain, 
        [τ, x], 
        [q(τ, x), c(τ, x)]
    )
    return pde_system
end

println("Start")
println()

pde_system = test(20.0, 60.0)    

x, τ = pde_system.ivs

q, c = pde_system.dvs

dx = 0.1

println("init discretizer")
println()
discretization = MOLFiniteDifference([x => dx], τ, approx_order = 2)

println("Start disc")
println()
prob = discretize(pde_system, discretization)

println("Start solve")
println()
# Solve ODE problem
sol = solve(prob, Tsit5(), saveat=0.2)

# Extract discrete values for x and τ
discrete_x = unique(sol[x])
discrete_τ = unique(sol[τ])

# Plot 3D surface for q(τ, x)
surface(discrete_x, discrete_τ, sol[q(τ, x)], xlabel="x", ylabel="τ", zlabel="q(τ,x)", title="q Concentration Profile")

# Plot 3D surface for c(τ, x)
surface(discrete_x, discrete_τ, sol[c(τ, x)], xlabel="x", ylabel="τ", zlabel="c(τ,x)", title="c Concentration Profile")
