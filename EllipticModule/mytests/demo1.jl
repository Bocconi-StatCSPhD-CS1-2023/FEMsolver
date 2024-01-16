using Plots

script_dir = @__DIR__
cd(script_dir)

include(joinpath(@__DIR__, "..", "src\\EllipticModule.jl"))

import .EllipitcPDE
const EP = EllipitcPDE

# The meshes available are from mesh0 to mesh4, the amplitude of the triangles that compose them change 
# (from the biggest to the smallest)

meshdir=joinpath(@__DIR__, "meshes\\mesh2")

#### EXAMPLE 1 ####

# Definition of the equation to solve. Terms not present don't need to be initated as they are missing by default.
Diff=1
f(x, y) = sin(π * x) * cos(y * π / 2) * π^2 * 5 / 4  # Forcing function

# boundary conditions
neuW(y) = π * cos(y * π / 2)  # Neumann in x = -1
neuE(y) = -π * cos(y * π / 2)  # Neumann in x = 1
dirS(x) = 0 # Dirichlet in y = 1
dirN(x) = 0 # Dirichlet in y = -1

# Numerical solution
coord, NumSol = EP.solution(meshdir;  Diff=Diff, f=f, neuW=neuW, neuE=neuE, dirN=dirN, dirS=dirS)
plot_numerical = surface(coord[:, 1], coord[:, 2], NumSol, c=:blues, title = "Numerical solution")

# Analytical solution
ureal(x, y) = sin(π * x) * cos(y * π / 2)  # This is the solution of my PDE
sol = ureal.(coord[:, 1], coord[:, 2]);
plot_analytical = surface(coord[:, 1], coord[:, 2], sol, seriestype = :surface, c=:reds, title = "Analytical solution")

# Comparison of the two results
combined_plot = plot(plot_numerical, plot_analytical, layout=(1, 2))
display(combined_plot)

# Computation of the L2 norm of the error 
error=EP.Err_L2_project(meshdir, NumSol, sol)
println("L2 Error for the given mesh is: $error")

