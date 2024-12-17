# 2D diffusion simulation in a square by using a degenerated form of the Allen-Cahn equation

## Statement of the problem
This example code solves the Allen-Cahn equation in a 2D domain $\Omega=[0,R]\times[0,R]$ with $R=1$ (inline mesh). 

 This test aims at assessing the use of an explicit scheme and the consistency of the diffusive part of the AllenCahn equation. 
 The Allen-Cahn equation is therefore simplified to degenerate into a diffusion equation and both equations are solved. The $L^2$ norm of the error is computed and must be exactly the same for both equations.



### __Governing equation__
Let us consider the following set of governing equations:

$$
\frac{\partial \phi}{\partial t}=M_\phi[\nabla \cdot{} \lambda \nabla \phi -\omega W'(\phi)]\text{ in }\Omega 
$$


$$
W(\phi)=\phi^2(1-\phi)^2
$$


The initial condition is simply a regularized Heavide centered at $R/2$.

$$\phi = \frac{1}{2}+\frac{1}{2}\tanh\bigg[2. \frac{(x - (R/2))}{ \epsilon_0}\bigg]$$

Particular solution for this problem are given in the following form:

$$\phi=\frac{1}{2} + \frac{1}{2} erf\bigg(\frac{x-L/2}{L_c} \bigg) $$ 

where $L_c=\sqrt{4M_\phi t}$


### __Boundary conditions__

Neumann boundary conditions are prescribed on $\Gamma_{lower}$ (y=0), $\Gamma_{upper}$ (y=R) and  $\Gamma_{right}$ (x=R):
$$
{\bf{n}} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma_{lower}  \cup \Gamma_{upper}  \cup \Gamma_{right}
$$


Dirichlet boundary conditions are prescribed on $\Gamma_{left}$ (x=0) :

$$
\phi=0 \text{ on }\Gamma_{left} 
$$



### __Parameters__
For this test, the following parameters are considered:

| Parameter                          | Symbol     | Value                       |
| ---------------------------------- | ---------- | --------------------------- |
| mobility coefficient               | $M_\phi$   | $1.e-2$                     |
| energy gradient coefficient        | $\lambda$  | $1$                         |
| surface tension                    | $\sigma$   | $1$                      |
| interface thickness                | $\epsilon$ | $1$                     |
| depth of the double-well potential | $\omega$   | $0$     |


### __Numerical scheme__

- Time marching: Euler Implicit scheme, $t\in[0,1]$, $\delta t=0.5(\delta x)^2/4M_\phi$
- Spatial discretization: uniform mesh with $\delta x=1/NN$ where $NN=20$ nodes in each direction
- Double-Well potential: implicit scheme (not important)


## Running

### __Using the binary__
```shell
./AllenCahn2Dtest7
```

### __Using ctest__

```shell
ctest -R AllenCahn2Dtest7
```


## References

None

## Intellectual Property

See [About page](../../../../../about.html) 