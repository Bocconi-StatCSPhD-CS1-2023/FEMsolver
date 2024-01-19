## FEM solution of Elliptic PDEs

The purpose of this project is to solve stationary PDEs of second order 
with $P_1$ Galerking Finite Element Method. The PDEs are of the kind
$$-\text{div}(\text{Diff}\nabla u)+\text{div}(\beta u)+cu=f, \quad \text { in } \Omega
$$
with Dirichlet or Neumann boundary conditions on 
$\Gamma_D\cup\Gamma_N=\partial\Omega$; in particular:
$$u=f_d \quad \text { on } \quad \Gamma_{D},\quad \text{Diff}\nabla u\cdot\nu=g \quad \text { on } \quad \Gamma_{N}.$$
$\text{Diff}$ is the diffusion coefficient, $\beta: \Omega \rightarrow \mathbb{R}^d$ is a vector field, and $c: \Omega \rightarrow \mathbb{R}$ is the reaction coefficient. From the physical point of view, this equation may represent the transport of a solute dissolved in a fluid that moves with the velocity field $\beta(x)$ and undergoes a linear chemical transformation with rate $c$.

See more details in the pdf file for the theory behind the numerical method.

The module EllipticModule allows to choose the presence or absence of the terms in the PDE, the constants $\text{Diff},\nabla\beta,c$ and the boundary conditions, that can be of the kind Neumann and homogeneous/non-homogeneous Dirichlet.
The meshes on which the code run are defined on the square [-1,1]x[-1,1]. To find the mesh vertices on the edges, the code takes advantage of the fact that they are exactly on the sides of the square.

### How to run the code
Defined the problem to solve, it is enough to define the parameters present in the equation:
$\text{Diff}\in\mathbb{R}, \nabla\beta\in\mathbb{R}^2,c\in\mathbb{R}$ and the function  $f:\mathbb{R}^2\rightarrow\mathbb{R}$.

Then it is necessary to define the boundary conditions. For each edge of the domain there must be a choice between Neumann or Dirichlet condition. So 4 functions $g:\mathbb{R}\rightarrow\mathbb{R}$ of the kind DirN/NeuN, DirE/NeuE, DirS/NeuS, DirW/NeuW have to be defined.

If you also want to apply the stabilizing method, you have to choose the $\delta$ weight.

To conclude, you have to choose the mesh on which run the method, and save the path to it in the variable "meshdir".

Finally run the function solution(meshdir; Diff=Diff, delta=delta, beta=beta, c=c, dirW=dirW, dirE=dirE, dirN=dirN, dirS=dirS) with the parameters chosen. Notice that they are keyword arguments, so if for example $\beta$ is not present in the equation, just do not define it.

The function solution provides the Nx2 coord matrix that indicates the points on the plane where the solution is computed and the vector NumSol, where the n-th entry correspond to the value of the solution at the point n.


In the demo file, there are two differential equations solved.

### Demo 1

The equation solved is 
$$-\Delta u=f, \quad \text { in } \Omega$$
where 
$$f(x,y)=\sin(πx)\cos(\frac{π}{2}y) π^2 \frac{5}{4}$$
with boundary conditions
$$u=0 \quad \text{ on } \quad \Gamma_{D}=\{(x, y):-1 \leq x \leq 1, y=\pm 1\}$$


$$\frac{\partial u}{\partial n}=q_{1}=-\pi \cos \left(\frac{\pi}{2} y\right) \quad \text { on } \quad \Gamma_{N 1}=\{(x, y): x=1,-1 \leq y \leq 1\}$$

$$\frac{\partial u}{\partial n}=q_{2}=\pi \cos \left(\frac{\pi}{2} y\right) \quad \text { on } \quad \Gamma_{N 2}=\{(x, y): x=-1,-1 \leq y \leq 1\}.$$

A qualitative comparison with the analytical solution $u$ is done, where $u$ is
$$u(x, y)=\sin (\pi x) \cos \left(\frac{\pi}{2} y\right).$$

Moreover as the exact solution is known, it has been computed the L2 error between the numerical and analytical solution. This operation iterated on refined meshes allows to see the goodness of the numerical method, as the convergence order of the method is 2.

### Demo 2
The equation solved is 
$$\text{div}(D\nabla u)-\text{div}(\beta u)=0  \text{ for (x,y)} \in \Omega$$
$$u(\text{(x,y)})=1, \text{ (x,y)} \in \Gamma_1$$
$$u(\text{ (x,y)})=0, \text{ (x,y)} \in \Gamma_2$$
where $\Gamma_1=\{-1\}\times[-1,1] \cup [-1,0.3]\times\{-1\}$ and $\Gamma_2=\{1\}\times[-1, 1] \cup [-1,1]\times\{1\}\cup [0.3,1]\times\{-1\}$ and $\beta$ is fixed and constant along all the domain $\Omega$: $\beta=(1,3)^T.$

We analyze the case when $D$ is small with respect to $\beta$, as in this situation the approximation error may be extremely large. We study the contribution of the Streamline upwind stabilization method (SUS) [see major details in the pdf file]
The SUS idea is to add an artificial diffusion coefficient along the flow lines of the $\beta$ convection field, since it is only along the velocity vector field that the upwind idea is most effective.

There is a qualitative comparison between the solution with and without the stabilization.

NOTE: For each demo you need to restart julia because the defined functions do not overwrite each other
 
