using Plots


script_dir = @__DIR__
cd(script_dir)

include(joinpath(@__DIR__, "..", "src\\EllipticModule.jl"))

import .EllipitcPDE
const EP = EllipitcPDE

meshdir=joinpath(@__DIR__, "meshes\\mesh2")

#### EXAMPLE 1 ####

f(x, y) = sin(π * x) * cos(y * π / 2) * π^2 * 5 / 4  # Forcing function

neuW(y) = π * cos(y * π / 2)  # Neumann in x = -1
neuE(y) = -π * cos(y * π / 2)  # Neumann in x = 1
dirS(x) = 0 # Dirichlet in y = 1
dirN(x) = 0 # Dirichlet in y = -1

coord, NumSol = EP.solution(meshdir; Diff=1, f=f, neuW=neuW, neuE=neuE, dirN=dirN, dirS=dirS)
plot_numerical = surface(coord[:, 1], coord[:, 2], NumSol, c=:blues, title = "Numerical solution")

ureal(x, y) = sin(π * x) * cos(y * π / 2)  # This is the solution of my PDE
sol = ureal.(coord[:, 1], coord[:, 2]);
# Plot Analytical solution
plot_analytical = surface(coord[:, 1], coord[:, 2], sol, seriestype = :surface, c=:reds, title = "Analytical solution")

combined_plot = plot(plot_numerical, plot_analytical, layout=(1, 2))
display(combined_plot)

error=EP.Err_L2_project(meshdir, NumSol, sol)
println("L2 Error for the given mesh is: $error")

ccc = EP.MassBalance(meshdir, NumSol; Diff=1, f=f, neuW=neuW, neuE=neuE, dirN=dirN, dirS=dirS)
