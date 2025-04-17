# 2D Allen-Cahn simulation  in piece of pellet


## Statement of the problem
This example code solves the Allen-Cahn equation in a 2D piece of pellet $\Omega$ (of radius $R=0.00465$ and angle $\theta=\pi/4$) (GMSH mesh) with a constant enthalpy of melting $\Delta H$. 
The pellet contains a circular inclusion of liquid.

This test aims at covering the use of attributes in the definition of the variables. 



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

The goal is to start with the given circular inclusion of liquid and follow its propagation through the pellet.

### __Boundary conditions__

Neumann boundary conditions are prescribed on $\Gamma$:
$$
{\bf{n}} \cdot{} \lambda \nabla \phi=0 \text{ on }\Gamma
$$


### __Parameters__
For this test, the following parameters are considered:

| Parameter                          | Symbol     | Value                       |
|------------------------------------|------------|-----------------------------|
| Enthalphy of melting               | $\Delta H$ | $-7.e3$                     |
| mobility coefficient               | $M_\phi$   | $1.e-5$                     |
| energy gradient coefficient        | $\lambda$  | $\frac{3}{2}\sigma\epsilon$ |
| surface tension                    | $\sigma$   | $0.06$                      |
| interface thickness                | $\epsilon$ | $5.e-4$                     |
| depth of the double-well potential | $\omega$   | $12{\sigma}/{\epsilon}$     |


### __Numerical scheme__

- Time marching: Euler Implicit scheme, $t\in[0,10]$, $\delta t=0.25$
- Spatial discretization: built from GMSH
- Double-Well potential: implicit scheme




## Running

### __Using the binary__
```shell
./AllenCahn2Dtest10
```

### __Using ctest__

```shell
ctest -R AllenCahn2Dtest10
```

### __In case of code coverage analysis__

```shell
make AllenCahn2Dtest10_coverage
```


## Post-processing

(to be written)

## Files & Dependencies

Source file : `main.cpp`
Mesh file : `camembert2D.msh` (build from `camembert.geo`)

## References

None

## Intellectual Property

See [About page](../../../../../about.html) 