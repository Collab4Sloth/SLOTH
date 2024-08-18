# 2D Allen-Cahn simulation in a square: steady quadratic manufactured solution

## Statement of the problem
This example code solves the Allen-Cahn equation in a 2D domain $\Omega=[0,1]\times[0,1]$ (inline mesh), with a spatial source term built in order to recover a steady quadratic solution:

$$
u({\bf x},t)=xy(x - 1)  (y - 1)
$$

### __Governing equation__
Let us consider the following set of governing equations:

$$
0=a\Delta u - H \partial_uf_{dw}(u)+S({\bf x},t)\text{ in }\Omega 
$$

$$
f_{dw}(u)=u^2(1-u)^2
$$

where $S({\bf x})=2 u({\bf x})  (1 - 2  u({\bf x}))  (1 -u({\bf x})) - [2  x  (x - 1) + 2  y  (y - 1)]$ 



### __Boundary conditions__




Dirichlet boundary conditions are prescribed on $\Gamma$ :

$$
\phi=0 \text{ on }\Gamma 
$$

### __Parameters__
For this test, all parameters are set to one.

### __Numerical scheme__

- Newton solver
- Spatial discretization: uniform mesh with $N=4$ nodes in each direction

  
## Running

### __Using the binary__
```shell
./AllenCahn2Dtest4
```

### __Using ctest__

```shell
ctest -R AllenCahn2Dtest4
```


## Post-processing

A convergence analysis was carried out and the results are presented below.
<div style="text-align:center">
<img title="2D MMS" src="../../../../../img/convergence_test4_2D.png" alt="" width="500"></p>
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
python3 convergence_study.py -l 1.
```

## References

None

## Intellectual Property

See [About page](../../../../../about.html) 