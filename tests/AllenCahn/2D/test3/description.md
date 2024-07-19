# 2D Allen-Cahn simulation in a square : steady linear manufactured solution

## Statement of the problem
This example code solves the Allen-Cahn equation in a 2D domain $\Omega=[0,1]\times[0,1]$ (inline mesh), with a spatiotemporal source term built in order to recover a steady linear solution:

$$
u({\bf x},t)=2x
$$

### __Governing equation__
Let us consider the following set of governing equations:

$$
\frac{\partial \phi}{\partial t}=\Delta \phi - W'(\phi)+S({\bf x},t)\text{ in }\Omega 
$$

$$
W(\phi)=\phi^2(1-\phi)^2
$$

where $S({\bf x},t)=2 u({\bf x},t)  (1 - 2u({\bf x},t)) (1 - u({\bf x},t))$ 



### __Boundary conditions__

Neumann boundary conditions are prescribed on $\Gamma_{lower}$ (y=0) and $\Gamma_{upper}$ (y=R):
$$
{\bf{n}} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma_{lower}  \cup \Gamma_{upper}
$$


Dirichlet boundary conditions are prescribed on $\Gamma_{left}$ :

$$
\phi=0 \text{ on }\Gamma _{left}
$$

Dirichlet boundary conditions are prescribed on $\Gamma_{right}$ :

$$
\phi=2 \text{ on }\Gamma _{right}
$$

### __Parameters__
For this test, all parameters are set to one.

### __Numerical scheme__

- Time marching: Euler Implicit scheme, $t\in[0,1]$, $\delta t=0.25$
- Spatial discretization: uniform mesh with $N=2$ nodes in each direction
- Double-Well potential: implicit scheme

  
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

(to be written)

## Files & Dependencies

Source file : `main.cpp`

## References

None

## Intellectual Property

See [About page](../../../../../about.html) 