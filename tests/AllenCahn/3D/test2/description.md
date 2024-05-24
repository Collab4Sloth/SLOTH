# 3D Allen-Cahn simulation  in piece of pellet


## Statement of the problem
This example code solves the Allen-Cahn equation in a 3D piece of pellet $\Omega$ (of radius $R=0.00465$, height $h=0.01$ and angle $\theta=\pi/4$) (GMSH mesh) with a constant enthalpy of melting $\Delta H$.

<div style="text-align:center">
<img  title="3D Geometry built in GMSH" src="../../../../../img/geometry_AC_3D_test2.png" alt="" width="500"></p> 
</div>

### __Governing equation__
Let us consider the following set of governing equations:
$$
\frac{\partial \phi}{\partial t}=M_\phi[\nabla \cdot{} \lambda \nabla \phi -\omega W'(\phi)+p(\phi)\Delta H]\text{ in }\Omega 
$$

$$
W(\phi)=\phi^2(1-\phi)^2
$$

$$
p(\phi)=\phi^3(6\phi^2-15\phi+10)
$$

The goal is to start with the given 3D hyperbolic tangent solution and follow its propagation through the pellet:

$$\phi({\bf{r}},t=t_{end}) = \frac{1}{2}+\frac{1}{2}\tanh\bigg[2. \frac{(x - (R/2))}{ \epsilon}\bigg]$$



### __Boundary conditions__

Neumann boundary conditions are prescribed on $\Gamma$:
$$
{\bf{n}} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma
$$


### __Parameters__
For this test, the following parameters are considered:

| Parameter                          | Symbol     | Value                       |
| ---------------------------------- | ---------- | --------------------------- |
| Enthalphy of melting               | $\Delta H$ | $7.e3$                      |
| mobility coefficient               | $M_\phi$   | $1.e-5$                     |
| energy gradient coefficient        | $\lambda$  | $\frac{3}{2}\sigma\epsilon$ |
| surface tension                    | $\sigma$   | $0.06$                      |
| interface thickness                | $\epsilon$ | $5.e-4$                     |
| depth of the double-well potential | $\omega$   | $12{\sigma}/{\epsilon}$     |


### __Numerical scheme__

- Time marching: Euler Implicit scheme, $t\in[0,0.25]$, $\delta t=0.25$ (one iteration only for covering test)
- Spatial discretization: built from GMSH
- Double-Well potential: implicit scheme

  

## Input file description
<iframe src="../../../../../../html/3D_2test2_2main_8cpp.html"  style="width: 100%; height: 100vh; border: none;"></iframe>


## Running

### __Using the binary__
```shell
./AllenCahn3Dtest2
```

### __Using ctest__

```shell
ctest -R AllenCahn3Dtest2
```

### __In case of code coverage analysis__

```shell
make AllenCahn3Dtest2_coverage
```


## Post-processing

(to be written)

## Files & Dependencies

Source file : `main.cpp`
Mesh file : `camembert3D.msh` (build from `camembert3D.geo`)

## References

None

## Intellectual Property

See [About page](../../../../../about.html) 
