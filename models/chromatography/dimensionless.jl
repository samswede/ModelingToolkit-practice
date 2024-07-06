using MethodOfLines, ModelingToolkit, DomainSets
using DifferentialEquations, Plots  

function test2()
    @parameters τ, ξ
    @variables q̃(..), c̃(..)  

    ∂τ = Differential(τ) 
    ∂ξ = Differential(ξ)
    ∂ξξ = Differential(ξ)^2
    
    ϵ = 0.4436 
    b = 0.0258
    k_f = 0.89
    D = 2.69e-4
    u = 0.0424
    qm = 76.85
    c0 = 8.8
    L = 20.0

    eq = [
        0 ~ - D / (L * u) * ∂ξξ(c̃(τ, ξ)) + ∂τ(c̃(τ, ξ)) + (1 - ϵ) / ϵ * qm / c0 * ∂τ(q̃(τ, ξ)) + ∂ξ(c̃(τ, ξ)),
        0 ~ u / L * ∂τ(q̃(τ, ξ)) + k_f * (q̃(τ, ξ) - b * c0 * c̃(τ, ξ) / (1 + b * c0 * c̃(τ, ξ)))
    ]

    domain = [ξ ∈ Interval(0.0, 1.0),
            τ ∈ Interval(0.0, 1.0)
    ]

    ic_bc = [
        # initial conditions
        q̃(0.0, ξ) ~ 0,
        c̃(0.0, ξ) ~ 0,
        # boundary conditions
        ∂ξ(c̃(τ, 0.0)) ~ u * L / D * (c̃(τ, 0.0) - 1),
        ∂ξ(c̃(τ, 1.0)) ~ 0.0
    ]

    @named pde_system = PDESystem(
        eq, ic_bc, domain, 
        [τ, ξ], 
        [q̃(τ, ξ), c̃(τ, ξ)]
    )
    return pde_system
end

pde_system2 = test2()    

τ, ξ = pde_system2.ivs

q̃, c̃ = pde_system2.dvs

dξ = 0.01

discretization = MOLFiniteDifference([ξ => dξ], τ, approx_order = 2)

println()
println("Created discretisation scheme")


prob2 = discretize(pde_system2, discretization)

println()
println("Discretised problem")


sol = solve(prob2, saveat = 1.0)

println()
println("Solved ODE problem")
