#import Pkg; Pkg.add("ModelingToolkitDesigner"); Pkg.add("ModelingToolkitStandardLibrary"); Pkg.add("GLMakie")
using ModelingToolkit
using ModelingToolkitDesigner
using GLMakie

using ModelingToolkitStandardLibrary
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC
import ModelingToolkitStandardLibrary.Blocks as B


@parameters t

@component function system(; name)

    pars = []

    systems = @named begin
        stp = B.Step(;height = 10e5, start_time = 0.005)
        src = IC.InputSource(;p_int=0)
        vol = IC.FixedVolume(;p_int=0, vol=10.0)
        res = IC.Pipe(5; p_int=0, area=0.01, length=500.0)
    end
   
    eqs = Equation[]
    
    ODESystem(eqs, t, [], pars; name, systems)
end


@named sys = system()

println("System created")


path = joinpath(@__DIR__, "design") # folder where visualization info is saved and retrieved
design = ODESystemDesign(sys, path);

println("Design created")

ModelingToolkitDesigner.view(design)
