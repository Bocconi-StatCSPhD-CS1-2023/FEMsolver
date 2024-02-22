using Plots

script_dir = @__DIR__
cd(script_dir)

include(joinpath(@__DIR__, "..", "src","EllipticModule.jl"))

import .EllipitcModule
const EP = EllipitcModule

# The meshes available are from mesh0 to mesh4, the amplitude of the triangles that compose them change 
# (from the biggest to the smallest)
meshdir=joinpath(@__DIR__, "meshes","mesh2")

#### EXAMPLE 2

Diff = 0.01
c = 0
beta = [1,3]

# boundary conditions
dirN(x) = 1 # Dirichlet in y = -1
dirE(x) = 0 # Dirichlet in y = 1
dirS(x) = 0
function dirW(x)
    return x > 0.3 ? 1 : 0
end

# Numerical solution without using the stabilization method UPU
delta = 0 
coord, NumSol = EP.solution(meshdir; Diff=Diff, delta=delta, beta=beta, c=c, dirW=dirW, dirE=dirE, dirN=dirN, dirS=dirS)
plot_numerical = surface(coord[:, 1], coord[:, 2], NumSol, c=:blues, title = "No stabilization")

# Numerical solution using the stabilization method UPU
delta = 0.1 
coord, NumSol = EP.solution(meshdir; Diff=Diff, delta=delta, beta=beta, c=c, dirW=dirW, dirE=dirE, dirN=dirN, dirS=dirS)
plot_numerical_st = surface(coord[:, 1], coord[:, 2], NumSol, c=:blues, title = "With stabilization")

# Comparison of the two results
combined_plot = plot(plot_numerical, plot_numerical_st, layout=(1, 2))
display(combined_plot)



