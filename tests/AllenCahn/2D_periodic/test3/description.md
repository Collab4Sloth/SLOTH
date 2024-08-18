# 2D Allen-Cahn simulation in a square: steady quadratic manufactured solution with with two periodic surface

## Statement of the problem
This example code solves the Allen-Cahn equation in a 2D domain $\Omega=[0,\pi]\times[0,\pi]$ (inline mesh), with a spatial source term built in order to recover a steady quadratic solution:

$$
u({\bf x},t)=\sin^2(2x) \cos^2(3y)
$$

### __Governing equation__
Let us consider the following set of governing equations:

$$
a\Delta u - H\partial_uf_{dw}(u)=-S({\bf x},t)\text{ in }\Omega 
$$

$$
f_{dw}(u)=u^2(1-u)^2
$$

where 

$$
\begin{multline}
S({\bf x})= 2 H \left(- 2 \sin^{2}{\left(2 x \right)} \cos^{2}{\left(3 y \right)} + 1\right) \left(- \sin^{2}{\left(2 x \right)} \cos^{2}{\left(3 y \right)} + 1\right) \sin^{2}{\left(2 x \right)} \cos^{2}{\left(3 y \right)} \\ - a \left(18 \sin^{2}{\left(2 x \right)} \sin^{2}{\left(3 y \right)} - 26 \sin^{2}{\left(2 x \right)} \cos^{2}{\left(3 y \right)} + 8 \cos^{2}{\left(2 x \right)} \cos^{2}{\left(3 y \right)}\right)
\end{multline}
$$


### __Boundary conditions__



Boundary boundary conditions are prescribed between $\Gamma_{lower}$ (y=0) and $\Gamma_{upper}$ (y=$\pi$).

Boundary boundary conditions are prescribed bteween  $\Gamma_{left}$ (x=0) and  $\Gamma_{right}$ (x=$\pi$).

### __Parameters__
For this test, all parameters are set to one.

### __Numerical scheme__

- Newton solver
- Spatial discretization: uniform mesh with $N=8$ nodes in each direction

  
## Running

### __Using the binary__
```shell
./AllenCahn2DPeriodictest3
```

### __Using ctest__

```shell
ctest -R AllenCahn2DPeriodictest3
```


## Post-processing

A convergence analysis was carried out and the results are presented below.
<div style="text-align:center">
<img title="2D MMS with periodic surface" src="../../../../../img/convergence_test3_2D_periodic.png" alt="" width="500"></p>
</div>

To carry out a convergence study, mesh parameters and element order must be defined as follows :

```c
std::vector<int> vect_NN{Nx_1,Nx_2,...,Nx_n};
std::vector<int> vect_order{o_1,o_2,...,o_m};
```

## Files & Dependencies

Source file : `main.cpp`

The convergence graph can be automatically generated using a python script such as :
```sh
python3 convergence_study.py -l 3.14159265359
```
## References

None

## Intellectual Property

See [About page](../../../../../about.html) 