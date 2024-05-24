# 3D Allen-Cahn simulation in a cube

## Statement of the problem
This example code solves the Allen-Cahn equation in a 3D domain $\Omega=[0,R]^3$ with $R=10^{-3}$ (inline mesh).

### __Governing equation__
Let us consider the following set of governing equations:

$$
\frac{\partial \phi}{\partial t}=M_\phi[\nabla \cdot{} \lambda \nabla \phi -\omega W'(\phi)]\text{ in }\Omega 
$$

$$
W(\phi)=\phi^2(1-\phi)^2
$$

The goal is to recover the 3D hyperbolic tangent solution

$$\phi({\bf{r}},t=t_{end}) = \frac{1}{2}+\frac{1}{2}\tanh\bigg[2. \frac{(x - (R/2))}{ \epsilon}\bigg]$$

with ${\bf{r}}=\sqrt{x^2+y^2+z^2}$, starting from hyperbolic tangent solution with a thinner interface $\epsilon_0=\epsilon/10$:

$$\phi({\bf{r}},t=0) = \frac{1}{2}+\frac{1}{2}\tanh\bigg[2. \frac{(x - (R/2))}{ \epsilon_0}\bigg]$$

### __Boundary conditions__

Neumann boundary conditions are prescribed on $\Gamma_{lower}$ (z=0), $\Gamma_{upper}$ (z=R), $\Gamma_{rear}$ (y=0) and $\Gamma_{front}$ (y=R):
$$
{\bf{n}} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma_{lower}  \cup \Gamma_{upper}  \cup 
\Gamma_{front}  \cup \Gamma_{rear} 
$$


Dirichlet boundary conditions are prescribed on $\Gamma_{left}$ (x=0) :

$$
\phi=0 \text{ on }\Gamma_{left} 
$$


Dirichlet boundary conditions are prescribed on $\Gamma_{right}$ (x=R):

$$
\phi=1 \text{ on }\Gamma_{right}
$$

### __Parameters__
For this test, the following parameters are considered:

| Parameter                          | Symbol     | Value                       |
| ---------------------------------- | ---------- | --------------------------- |
| mobility coefficient               | $M_\phi$   | $1.e-5$                     |
| energy gradient coefficient        | $\lambda$  | $\frac{3}{2}\sigma\epsilon$ |
| surface tension                    | $\sigma$   | $0.06$                      |
| interface thickness                | $\epsilon$ | $5.e-4$                     |
| depth of the double-well potential | $\omega$   | $12{\sigma}/{\epsilon}$     |


### __Numerical scheme__

- Time marching: Euler Implicit scheme, $t\in[0,0.25]$, $\delta t=0.25$ (one iteration only for covering test)
- Spatial discretization: uniform mesh with $N=30$ nodes in each direction
- Double-Well potential: implicit scheme

  

## Input file description
<iframe src="../../../../../../html/3D_2test1_2main_8cpp.html"  style="width: 100%; height: 100vh; border: none;"></iframe>


## Running

### __Using the binary__
```shell
./AllenCahn3Dtest1
```

### __Using ctest__

```shell
ctest -R AllenCahn3Dtest1
```

### __In case of code coverage analysis__

```shell
make AllenCahn3Dtest1_coverage
```


## Post-processing

(to be written)

## Files & Dependencies

Source file : `main.cpp`

## References

None

## Intellectual Property

See [About page](../../../../../about.html) 