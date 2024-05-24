# 2D Allen-Cahn simulation in a periodic domain

## Statement of the problem
This example code solves the Allen-Cahn equation in a periodic 2D domain $\Omega=[0,R]\times[0,R]$ with $R=2\pi$ (inline mesh). 

The periodicity is here applied after building the mesh from scratch using MFEM functionalities but a periodic mesh built with GMSH can be used.

This case is based on :

```
GUO, Ruihan, JI, Liangyue, et XU, Yan. High order local discontinuous Galerkin methods for the Allen-Cahn equation: analysis and simulation. Journal of Computational Mathematics, 2016, p. 135-158.
```

### __Governing equation__
Let us consider the following set of governing equations:

$$
\frac{\partial \phi}{\partial t}=M_\phi[\nabla \cdot{} \lambda \nabla \phi -\omega W'(\phi)] + S \text{ in }\Omega 
$$

$$
W(\phi)=\phi^3 - \phi
$$


The goal is to recover the sinusoide solution

$$\phi({\bf{r}},t=t_{end}) = e^{-2. t_{end}} sin(x+y)$$

by prescribing the following source term

$$ S({\bf{r}},t) = u({\bf{r}},t)^3 - u({\bf{r}},t)$$

$$u({\bf{r}},t)= e^{-2. t} sin(x+y)$$ls

### __Boundary conditions__

Boundary boundary conditions are prescribed between $\Gamma_{lower}$ (y=0) and $\Gamma_{upper}$ (y=R).

Boundary boundary conditions are prescribed bteween  $\Gamma_{left}$ (x=0) and  $\Gamma_{right}$ (x=R).

### __Parameters__
For this test, the following parameters are considered:

| Parameter                          | Symbol     | Value               |
| ---------------------------------- | ---------- | ------------------- |
| mobility coefficient               | $M_\phi$   | $5.e-5$             |
| surface tension                    | $\sigma$   | $0.06$              |
| energy gradient coefficient        | $\lambda$  | $3\sigma\epsilon/2$ |
| interface thickness                | $\epsilon$ | $3.e-4$             |
| depth of the double-well potential | $\omega$   | $12\sigma/\epsilon$ |


### __Numerical scheme__

- Time marching: Euler Implicit and Explicit  schemes, $t\in[0,0.2]$, $\delta t=0.05$
- Spatial discretization: uniform mesh with $N=32$ nodes in each direction
- Double-Well potential: implicit scheme

  

## Input file description
<iframe src="../../../../../../html/2D__periodic_2test1_2main_8cpp.html"  style="width: 100%; height: 100vh; border: none;"></iframe>


## Running

### __Using the binary__
```shell
./AllenCahn2D_periodictest2
```

### __Using ctest__

```shell
ctest -R AllenCahn2D_periodictest2
```

### __In case of code coverage analysis__

```shell
make AllenCahn2D_periodictest_2_coverage
```


## Post-processing

(to be written)

## Files & Dependencies

Source file : `main.cpp`

## References

None

## Intellectual Property

See [About page](../../../../../about.html) 