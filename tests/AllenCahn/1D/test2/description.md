# 1D Allen-Cahn simulation  along a radius with constant melting factor


## Statement of the problem
This example code solve one Allen-Cahn equations in a 1D domain $\Omega=[0,L]$ with $L=10^{-3}$ (inline mesh) with a constant melting factor (ie. constant driving force). 


### __Governing equation__

  Let us consider the following set of governing equations:

$$
\frac{\partial \phi}{\partial t}=M_\phi[\nabla \cdot{} \lambda \nabla \phi -\omega W'(\phi) + p'(\phi)\mathcal{R}]\text{ in }\Omega 
$$

$$
W(\phi)=\phi^2(1-\phi)^2
$$

$$
p(\phi) = \phi^3(6\phi^2-15\phi+10)
$$


starting from hyperbolic tangent solution :

$$\phi(r,t=0) = \frac{1}{2}+\frac{1}{2}\tanh\bigg[2. \frac{(r - L/2)}{ \epsilon}\bigg]$$

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
| melting factor                     | $\mathcal{R}$| $7.e-3$                   |

### __Numerical scheme__

- Time marching: Euler Implicit/Euler Explicit/ RK4 schemes, $t\in[0,10]$, $\delta t=0.1$
- Spatial discretization: uniform mesh with $N=30$ nodes 
- Double-Well potential: implicit scheme
## Input file description
<iframe src="../../../../../../html/1D_2test1_2main_8cpp.html"  style="width: 100%; height: 100vh; border: none;"></iframe>


## Running

### __Using the binary__
```shell
./AllenCahn1Dtest1
```

### __Using ctest__

```shell
ctest -R AllenCahn1Dtest2
```

### __In case of code coverage analysis__

```shell
make AllenCahn1Dtest2_coverage
```


## Post-processing

This test highlights the ability of post-processing to extract the position of an isovalue from the solution. In this case, the position $\phi$=0.5 is extracted :
<div style="text-align:center">
<img title="Interface position for 1D simulation" src="../../../../../img/AC_1D_test2.png" alt="" width="500"></p>
</div>

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
