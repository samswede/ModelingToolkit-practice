# The purpose of this model is to implement the method of lines for the convection-diffusion equation without the MOL library.
# Instead, we will use ModelingToolkit.jl to symbolically derive the spatial discretization and then use DifferentialEquations.jl to solve the resulting ODE system.

