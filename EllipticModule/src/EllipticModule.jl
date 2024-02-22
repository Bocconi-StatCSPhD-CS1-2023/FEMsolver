module EllipitcModule

using DelimitedFiles
using LinearAlgebra
using SparseArrays

function input_data(meshdir)
    # Import the mesh that triangularizes the domain and says the boundary nodes. The domain is the square [-1,1]x[-1,1]
    # Each vertex is enumerated, triang is a Nelemx3 matrix where the i^th row corresponds to the vertices that define the i^th element
    # Coord is a Nodesx2 matrix where the i^th row corresponds to the coordinates of the i^th node

    # meshdir is the path of the directory where the meshes are saved
    meshtriang = string(joinpath(meshdir, "triang.dat"))
    meshcoord = string(joinpath(meshdir, "xy.dat"))

    triang = readdlm(meshtriang) #<CB># you can read integers adding `Int` as the second argument
    xy = readdlm(meshcoord)

    triang = convert(Matrix{Int64},triang[:, 1:3])
    coord = xy[:, 1:2]

    Nelem = size(triang, 1)
    Nodes = size(coord, 1)

    boundary_nodes = Dict("N" => Vector{Int}(), "S" => Vector{Int}(), "E" => Vector{Int}(), "W" => Vector{Int}())
    for j in 1:size(coord,1)
        if coord[j, 2] == -1
            push!(boundary_nodes["S"], j)
        end
        if coord[j, 2] == 1
            push!(boundary_nodes["N"], j) 
        end
        if coord[j, 1] == -1
            push!(boundary_nodes["W"], j)
        end
            if coord[j, 1] == 1
            push!(boundary_nodes["E"], j)
        end
    end

    return Nelem, Nodes, triang, coord, boundary_nodes
end


function localBasis(Nelem, triang, coord)
    # returns the Bloc and Cloc Nelem x 3 matrices where in the i^th row there are the coefficients b_1,b_2b_3 of the linear basis functions 
    # that have support in the element
    # Area is a Nelem vector where there are listed the areas of the triangles 
    Area = zeros(Nelem)
    Bloc = zeros(Nelem, 3)
    Cloc = zeros(Nelem, 3)

    B = ones(3, 3)

    for i in 1:Nelem
        for j in 1:3
            for k in 1:2
                B[j, k + 1] = coord[triang[i, j], k]
            end
        end

        for m in 1:3
            tnoto = zeros(3)
            tnoto[m] = 1
            sol = B \ tnoto
            Bloc[i, m] = sol[2]
            Cloc[i, m] = sol[3]
        end

        Area[i] = 0.5 * abs(det(B))
    end

    return Bloc, Cloc, Area
end


function localLengthBoundary(Nelem, triang, coord, Area, boundary_nodes)
    # The function returns 
    # ElementBoundary: dictionary where the elements that have vertices on the boundary are listed
    # Length: a Nelem vector with listed the length of the eventual boundary edge

    ElementBoundary = Dict("N" => Vector{Int}(), "S" => Vector{Int}(), "E" => Vector{Int}(), "W" => Vector{Int}())
    Length = zeros(Nelem)

    for i in 1:Nelem

        for key in keys(boundary_nodes)

            inters = intersect(triang[i,:], boundary_nodes[key])
            if length(inters) == 2
                push!(ElementBoundary[key], i)

                if coord[inters[1],1]==coord[inters[2],1]
                    ll = abs( coord[inters[1],2]-coord[inters[2],2])
                else 
                    ll = abs( coord[inters[1],1]-coord[inters[2],1])
                end

                Length[i] = ll
            end
        end
    end

    return ElementBoundary, Length
end


function stiffBuild(Nelem, Nodes, triang, Bloc, Cloc, Area, Diff)
    # Assembly of the stiffness matrix from local contributions (2D P1 Galerkin)
    # stiff_{ij}=\int Diff \phi'_i\phi'_j
    # Diff=diffusion coefficient
    # OUTPUT: stiffMat matrix Nelem x Nelem
    # The generic element a_ij=int_T a_i*a_j+b_i*b_j=area(T)[a_i*a_j+b_i*b_j]

    stiffMat = sparse(zeros(Nodes, Nodes)) #<CB># Use spzeros instead of this!
    for iel in 1:Nelem
        area = Area[iel]
        for iloc in 1:3
            iglob = triang[iel, iloc]
            for jloc in 1:3
                jglob = triang[iel, jloc]
                # for every node there are more contributions to sum
                stiffMat[iglob, jglob] += Diff * area * (Bloc[iel, iloc] * Bloc[iel, jloc] + Cloc[iel, iloc] * Cloc[iel, jloc])
            end
        end
    end

    return stiffMat
end

function Transport_build(Nelem, Nodes, triang, Bloc, Cloc, Area, beta)
    # Assembly of the matrix B from local contributions (2D P1 Galerkin)
    # B_{ij} = \int ( beta \cdot \phi'_j ) \phi_i
    # The generic element B_ij = \int_T ( (beta_1 b_j) + (beta_2 c_j))(a_i + b_i x + c_i y)
    # Using the trapezoidal rule: B_ij = \int_T ((beta_1 b_j) + (beta_2 c_j))(a_i + b_i x + c_i y) = (beta_1 b_j + beta_2 c_j)|T|/3

    B = sparse(zeros(Nodes, Nodes)) # a lot of zeros, convenient

    for iel in 1:Nelem
        area = Area[iel] / 3
        for iloc in 1:3
            iglob = triang[iel, iloc]
            for jloc in 1:3
                jglob = triang[iel, jloc]
                B[iglob, jglob] += (Bloc[iel, jloc] * beta[1] + Cloc[iel, jloc] * beta[2]) * area
            end
        end
    end

    return B
end


function MassBuild(Nelem, Nodes, triang, Area, c)

    #Assembly of the mass matrix from local contributions (2D P1 Galerkin)
    # mass_{ij}=\int c \phi_i\phi_j
    #The generic element m_ij = area(T)[1/6 if i=j, 1/12 otherwise]

    MassMat = sparse(zeros(Nodes, Nodes))
    for iel in 1:Nelem
        area = Area[iel]
        for iloc in 1:3
            iglob = triang[iel, iloc]
            for jloc in 1:3
                jglob = triang[iel, jloc]
                # for every node there are more contributions to sum
                if iglob == jglob
                    MassMat[iglob, jglob] += c * area / 6
                else MassMat[iglob, jglob] += c * area / 12
                end
            end
        end
    end
    return MassMat
end

function Forcing_vector(Nelem, Nodes, coord, triang, Bloc, Cloc, Area, ffunc)
    # Using the TRAPEZOIDAL RULE for evaluating the integral
    # (it is exact on linear functions)
    # has a convergence error compatible with linear P1 FEM

    # FORCING COMPONENT
    rhs = zeros(Nodes)

    for iel in 1:Nelem
        area = Area[iel]
        for iloc in 1:3
            iglob = triang[iel, iloc]
            # have to consider the contribution to the node of all the elements
            # that have it as a vertex
            rhs[iglob] += ffunc(coord[iglob, 1], coord[iglob, 2]) * area / 3
        end
    end

    return rhs
end


function stabilization_build(Nelem, Nodes, triang, Bloc, Cloc, Area, beta, delta)
    # Initialization of the stabilization matrix AS
    AS = sparse(zeros(Nodes, Nodes)) 

    for iel in 1:Nelem
        area = Area[iel]
        for iloc in 1:3
            iglob = triang[iel, iloc]
            for jloc in 1:3
                jglob = triang[iel, jloc]
                AS[iglob, jglob] += (Bloc[iel, jloc] * beta[1] + Cloc[iel, jloc] * beta[2]) *
                                    (Bloc[iel, iloc] * beta[1] + Cloc[iel, iloc] * beta[2]) * area
            end
        end
    end
    
    return delta*AS
end


### BOUNDARY CONDITIONS

function forcing_Neumann(ElementBoundary, Length, boundary_nodes, coord, Nodes, triang, neuW, neuE, neuS, neuN)
    # NEUMANN COMPONENT
    rhs = zeros(Nodes)

    ##<CB># Lots of code duplication here! Instead of copy-pasting, do another loop over "W", "E", "S", "N"
    ##      Things that need changing within the loop can be stored in a dict; or perhaps you can do
    ##          for (neu, dir, c) in [(neuW, "W", -1), (newE, "E", 1), ...]
    if neuW != nothing
        NeuNod = boundary_nodes["W"]
        for iel in ElementBoundary["W"]
            area = Length[iel]
            for iloc in 1:3
                iglob = triang[iel, iloc]
                if coord[iglob, 1] == -1
                    rhs[iglob] += neuW(coord[iglob, 2]) * area / 2
                end
            end
        end
    end 
    if neuE != nothing
        NeuNod = boundary_nodes["E"]
        for iel in ElementBoundary["E"]
            area = Length[iel]
            for iloc in 1:3
                iglob = triang[iel, iloc]
                if coord[iglob, 1] == 1
                    rhs[iglob] += neuE(coord[iglob, 2]) * area / 2
                end
            end
        end
    end
    if neuS != nothing
        NeuNod = boundary_nodes["S"]
        for iel in ElementBoundary["S"]
            area = Length[iel]
            for iloc in 1:3
                iglob = triang[iel, iloc]
                if coord[iglob, 2] == -1
                    rhs[iglob] += neuS(coord[iglob, 1]) * area / 2
                end
            end
        end
    end
    if neuN != nothing
        NeuNod = boundary_nodes["N"]
        for iel in ElementBoundary["N"]
            area = Length[iel]
            for iloc in 1:3
                iglob = triang[iel, iloc]
                if coord[iglob, 2] == 1
                    rhs[iglob] += neuN(coord[iglob, 1]) * area / 2
                end
            end
        end
    end

    return rhs
end


function forcing_Dir(coord, Nodes, boundary_nodes, stiffMat, rhs, dirW, dirE,  dirN, dirS)
    # Impose Dirichlet BCs with the penalty method
    
    penalty = 1e15
    AA = copy(stiffMat)
    rh = copy(rhs)

    ##<CB># see above
    if dirW != nothing
        DirNod = boundary_nodes["W"]
        for iglob in DirNod
            AA[iglob, iglob] = penalty
            y = coord[iglob,2]
            rh[iglob] = penalty * dirW(y)
        end
    end
    if dirE != nothing
        DirNod = boundary_nodes["E"]
        for iglob in DirNod
            AA[iglob, iglob] = penalty
            y = coord[iglob,2]
            rh[iglob] = penalty * dirE(y)
        end
    end
    if dirN != nothing
        DirNod = boundary_nodes["N"]
        for iglob in DirNod
            AA[iglob, iglob] = penalty
            y = coord[iglob,1]
            rh[iglob] = penalty * dirN(y)
        end
    end
    if dirS != nothing
        DirNod = boundary_nodes["S"]
        for iglob in DirNod
            AA[iglob, iglob] = penalty
            y = coord[iglob,1]
            rh[iglob] = penalty * dirS(y)
        end
    end

    return AA, rh
end

function solution(meshdir; Diff= 0, delta=0, beta= zeros(2), c=0, f = (x,y)->0, neuW =nothing, neuE =nothing, neuS =nothing, neuN =nothing, 
                                                                            dirN =nothing, dirS =nothing, dirW =nothing, dirE =nothing  )

    Nelem, Nodes, triang, coord, boundary_nodes = input_data(meshdir)
    Bloc, Cloc, Area = localBasis(Nelem, triang, coord)
    ElementBoundary, Length = localLengthBoundary(Nelem, triang, coord, Area, boundary_nodes)

    S = stiffBuild(Nelem, Nodes, triang, Bloc, Cloc, Area, Diff)
    T = Transport_build(Nelem, Nodes, triang, Bloc, Cloc, Area, beta)
    M = MassBuild(Nelem, Nodes, triang, Area, c)
    L = stabilization_build(Nelem, Nodes, triang, Bloc, Cloc, Area, beta,delta)
    A = S + T + M + L

    rhs0 = Forcing_vector(Nelem, Nodes, coord, triang, Bloc, Cloc, Area, f)
    rhs1 = forcing_Neumann(ElementBoundary, Length, boundary_nodes, coord, Nodes, triang, neuW, neuE, neuS, neuN)
    rhs = rhs0 + rhs1
    AA, rh = forcing_Dir(coord, Nodes, boundary_nodes, A, rhs, dirW, dirE, dirN, dirS)

    NumSol = AA \ rh

    return coord, NumSol
end

function Err_L2_project(meshdir, uh, ureal)

    # using the trapezoidal rule for evaluating the integral
    # (it is exact on linear functions)
    # has a convergence error compatible with linear P1 FEM

    Nelem, Nodes, triang, coord, boundary_nodes = input_data(meshdir)
    Bloc, Cloc, Area = localBasis(Nelem, triang, coord)

    # Initialize error
    err = 0.0
    
    # Loop over elements
    for iel in 1:Nelem
        area = Area[iel]
        
        # Loop over local nodes
        for iloc in 1:3
            iglob = triang[iel, iloc]
            
            # Compute the squared difference and accumulate
            err += (uh[iglob] - ureal[iglob])^2 * area / 3
        end
    end
    
    # Take square root of accumulated error
    err = sqrt(err)
    
    return err
end


end #module EllipitcPDE
