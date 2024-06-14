# 1D Allen-Cahn simulation  along a radius


## Statement of the problem
This example code solves three  Allen-Cahn equations in a 1D domain $\Omega=[0,R]$ with $R=10^{-3}$ (inline mesh). 
Each equation corresponds to a problem consider inside the same coupling. There is no interaction between the equations (just for test).
Each equation uses a different numerical scheme.

### __Governing equation__
Let us consider the following set of governing equations:

$$
\frac{\partial \phi}{\partial t}=M_\phi[\nabla \cdot{} \lambda \nabla \phi -\omega W'(\phi)]\text{ in }\Omega 
$$

$$
W(\phi)=\phi^2(1-\phi)^2
$$

The goal is to recover the 1D hyperbolic tangent solution

$$\phi(r,t=t_{end}) = \frac{1}{2}+\frac{1}{2}\tanh\bigg[2. \frac{(r - (R/2))}{ \epsilon}\bigg]$$

starting from hyperbolic tangent solution with a thinner interface $\epsilon_0=\epsilon/10$:

$$\phi(r,t=0) = \frac{1}{2}+\frac{1}{2}\tanh\bigg[2. \frac{(r - (R/2)))}{ \epsilon_0}\bigg]$$

In this case, the initial condition is defined as a `std::function`.

### __Boundary conditions__

Neumann boundary conditions are prescribed on $\Gamma_{left}$ (r=0) and $\Gamma_{right}$ (r=R):
$$
{\bf{n}} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma_{left} \cup \Gamma_{right}
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

- Time marching: Euler Implicit/Euler Explicit/ RK4 schemes, $t\in[0,50]$, $\delta t=0.25$
- Spatial discretization: uniform mesh with $N=30$ nodes 
- Double-Well potential: implicit/explicit/SemiImplicit scheme

  

## Input file description
<iframe src="../../../../../../html/1D_2test1_2main_8cpp.html"  style="width: 100%; height: 100vh; border: none;"></iframe>


## Running

### __Using the binary__
```shell
./AllenCahn1Dtest1
```

### __Using ctest__

```shell
ctest -R AllenCahn1Dtest1
```

### __In case of code coverage analysis__

```shell
make AllenCahn1Dtest1_coverage
```


## Post-processing

(to be written)

## Files & Dependencies

Source file : `main.cpp`

## References

None

## Intellectual Property

See [About page](../../../../../about.html) 

<!-- !!! note
    Example of a note.
!!! tip "Custom title"
    Example tip. -->
