using MethodOfLines, ModelingToolkit, DomainSets, DifferentialEquations

function test(L::Float64, T::Float64)
    @parameters t, x
    @variables q(..), c(..)  

    ∂t = Differential(t) 
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
        0 ~ D * ∂xx(c(t, x)) + ∂t(c(t, x)) + (1 - ϵ) / ϵ * ∂t(q(t, x)) + u * ∂x(c(t, x)),
        0 ~ ∂t(q(t, x)) + k_f * (q(t, x) + qm * b * c(t, x) / (1 + b * c(t, x)))
    ]

    domain = [x ∈ Interval(0.0, L),
            t ∈ Interval(0.0, T)
    ]

    ic_bc = [
        # initial conditions
        q(0.0, x) ~ 0.0,
        c(0.0, x) ~ 0.0,
        # boundary conditions
        ∂x(c(t, 0)) ~ u * (c(t, 0) - c0) / D, #AH! c(t, x) should be c(t, 0) (I think)... he defined C_in as c0.
        c(t, L) ~ c(t, L - 0.1) # Approximating a Neumann BC, ∂x(c(t, L)) ~ 0
    ]

    @named pde_system = PDESystem(
        eq, ic_bc, domain, 
        [t, x], 
        [q(t, x), c(t, x)]
    )
    return pde_system
end

println("Start")
println()

pde_system = test(20.0, 60.0)    

println("Built pde system")
println()

x, t = pde_system.ivs

q, c = pde_system.dvs

dx = 0.1

discretization = MOLFiniteDifference([x => dx], t, approx_order = 2)

prob = discretize(pde_system, discretization)

println("Successfully discretized")
println()