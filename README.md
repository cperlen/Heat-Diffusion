# Heat-Diffusion

We consider the two-dimensional heat diffusion equation given by $$ \frac{\partial  T}{\partial t} = \kappa \nabla^2 T$$ with $\kappa = 1$ on the domain $0 \leq x \leq \pi,$ $0 \leq y \leq \pi$ subject to the boundary conditions \begin{align}T(x,0) &= \text{cos}^2(x) \\ T(x,\pi) &= \text{sin}^2(x) \\ T(0,y) &= T(\pi,y) \text{ (periodic in x)} \end{align}

We solve this through centered finite differences in space and the forward Euler method in time.  We implement this via

heat\_serial: Serial version.  Run with ./heat_serial <nx> for a solution with grid size $nx^2$

heat\_omp: Parallel version utilizing open mp.  Run with ./heat_omp <nx> <nthreads> for a solution with grid size $nx^2$ utilizing nthreads threads.  

• heat\_MPI: Parallel version utilizing MPI. Run with ./heat_mpi {nx}. 
