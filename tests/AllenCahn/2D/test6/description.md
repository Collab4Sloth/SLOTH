# 2D Allen-Cahn simulation : the PFHUB Benchmark

TO DO


## Statement of the problem
Test case from [PFHUB](https://pages.nist.gov/pfhub/benchmarks/benchmark7.ipynb/)


### __Governing equation__
Let us consider the following set of governing equations:
$$
\frac{\partial \phi}{\partial t}=\kappa \nabla^2 \phi -\omega f_{dw}'(\phi) + S_{MMS} \quad\text{   in }\Omega 
$$

$$
f_{dw}(\phi)=\phi^2(1-\phi)^2
$$




$$\phi_{\text{exact}}({\bf{r}},t=t_{end}) = \frac{1}{2}+\frac{1}{2}\tanh\bigg[\frac{y-\frac{1}{4}-\alpha_1t\sin{(\alpha_2x)}-\alpha_3\sin(\alpha_4x+\alpha_4t)}{\sqrt{2\kappa}}\bigg]$$



### __Boundary conditions__

Boundary conditions are prescribed between $\Gamma_{left}$ (y=0) and $\Gamma_{right}$ (y=$1$).
Neumann boundary conditions are prescribed on top and bottom:
$$
{\bf{n}} \cdot{} \nabla \phi=0 \text{ on }\Gamma
$$



### __Parameters__
cf. [PFHUB](https://pages.nist.gov/pfhub/benchmarks/benchmark7.ipynb/)


### __Numerical scheme__

- Time marching: Euler Implicit scheme, $t\in[0,0.1]$, $\delta t=0.01$
- Spatial discretization: quadrilater
- Double-Well potential: implicit scheme

  

## Input file description



## Running

### __Using the binary__
```shell
./AllenCahn2Dtest6
```

### __Using ctest__

```shell
ctest -R AllenCahn2Dtest6
```

### __In case of code coverage analysis__

```shell
make AllenCahn2Dtest6_coverage
```


## Files & Dependencies

Source file : `main.cpp`

## References

[ PFHub: The Phase Field Community Hub ](https://pages.nist.gov/pfhub/) 


## Intellectual Property

See [About page](../../../../../about.html) 