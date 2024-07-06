using Plots
function main()
    # Physical parameters
    L   = 20       # Length of bed ,cm
    ϵ   = 0.4436   # bed porosity
    ci  = 8.8      # Initial concentration ,mg/ml
    qm  = 76.85    # maximum adsorption capacity ,mg/ml
    b   = 0.0258   # Langmuir isotherm constant 
    u   = 0.0424   # Velocity ,cm/s
    k_f = 0.89     # mass transfer coefficient ,s-1  
    D   = 2.69e-4  # Axial diffusion coefficient ,cm2/s 
    # Numerical parameters
    ncx = 100
    Δx  = L/(ncx)
    Δt  = 3e-1
    nt  = 2000
    nit = 3
    # Grid and arrays allocations
    xc  = LinRange(Δx/2, L-Δx/2, ncx)
    c   = zeros(ncx+2) # includes 2 ghosts cells
    c0  = zeros(ncx+2)
    q   = zeros(ncx)
    q0  = zeros(ncx)
    # Time loop
    @views for it=1:nt
        @. c0 = c
        @. q0 = q
        # Iteration loop
        for iter = 1:nit
            # Boundaries
            c[end] = c[end-1]
            c[1]   = c[2] - Δx*u/D*(0.5*(c[2]) - ci)
            # Update balance laws
            @. c[2:end-1] = c0[2:end-1] + Δt*(D*(c[1:end-2] + c[3:end-0] - 2.0.*c[2:end-1])/Δx^2 - (1.0-ϵ)/ϵ*(q-q0)/Δt - (u>0)*u*(c[2:end-1] - c[1:end-2])/Δx)
            @. q          = q0          + Δt*(-k_f*(q - (qm*b*c[2:end-1])/(1.0 + b*c[2:end-1])) )
        end
        # Visualise
        if mod(it, 100) == 0
            p = plot(xlabel = "x", ylabel="c, q", ylims=(0,25))
            p = plot!(xc, c[2:end-1], label="c")
            p = plot!(xc, q, label="q")
            display(p)
            sleep(0.05)
        end
    end
end

main()