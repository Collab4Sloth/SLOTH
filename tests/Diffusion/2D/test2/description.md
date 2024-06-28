# 2D Heat transfer in a star (non linear case)


## Statement of the problem
This example code is based on the example 16 on MFEM's website. It solves a heat transfer equation in a 2D star $\Omega$ (GMSH mesh coming from MFEM examples).

<div style="text-align:center">
<img  title="2D Geometry built in GMSH" src="../../../../../img/geom_Diff_2D_test1.png" alt="" width="500"></p> 
</div>
 
### __Governing equation__
Let us consider the following set of governing equations:
$$
\frac{\partial T}{\partial t}=[\nabla \cdot{} \displaystyle(\kappa+\alpha T) \nabla T ]\text{ in }\Omega 
$$


The goal is to start with a fictitious temperature value set to $2$ inside a circle of radius 0.5, $1$ outside and follow heat diffusion through the star.

The diffusion coefficient implicitly depends on temperature leading to a non linear problem.

### __Boundary conditions__

Neumann boundary conditions are prescribed on $\Gamma$:
$$
{\bf{n}} \cdot{} \lambda \nabla T=0 \text{ on }\Gamma
$$


### __Parameters__
For this test, the following parameters are considered:

| Parameter             | Symbol   | Value   |
| --------------------- | -------- | ------- |
| Primary coefficient   | $\kappa$ | $0.5$   |
| Secondary coefficient | $\alpha$ | $1.e-2$ |


### __Numerical scheme__

- Time marching: Euler Implicit scheme, $t\in[0,0.5]$, $\delta t=0.1$
- Spatial discretization: built from GMSH + quadratic FE + 2 uniform levels of refinement 
- Double-Well potential: implicit scheme

<div style="text-align:center">
<img title="2D Finite Element Mesh built in GMSH" src="../../../../../img/mesh_Diff_2D_test1.png" alt="" width="500"></p>
</div>
 



## Running

### __Using the binary__
```shell
./Diffusion2Dtest1
```

### __Using ctest__

```shell
ctest -R Diffusion2Dtest2
```

### __In case of code coverage analysis__

```shell
make Diffusion2Dtest2_coverage
```


## Post-processing

(to be written)

## Files & Dependencies

Source file : `main.cpp`
Mesh file : `star2D.msh` (taken from MFEM examples)

## References

None

## Intellectual Property

See [About page](../../../../../about.html) 