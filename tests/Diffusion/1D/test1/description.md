# 1D semi-infinite domain simulation


## Statement of the problem
This example code is based on the analytical solution for a semi-infinite domain.


 
### __Governing equation__
Let us consider the following set of governing equations:
$$
\frac{\partial c}{\partial t}=[\nabla \cdot{} M\nabla c]\text{ in }\Omega 
$$



### __Boundary conditions__
Dirichlet boundary conditions are prescribed on $\Gamma_{left}$:

$$
c(x=0,t) = 1 
$$

Neumann boundary conditions are prescribed on $\Gamma_{right}$:

$$
{\bf{n}} \cdot{} M \nabla c=0 \text{ on }\Gamma_{right}
$$

### __Initial conditions__

The composition is assumed to be constant at $t=0$: 
$$
c(x,t=0) = 0 
$$

### __Analytical solution__

The analytical solution can be written:
$$
c(x,t) = 1-\text{erf}\left(\dfrac{x}{\sqrt{4Mt}}\right)
$$

### __Parameters__
For this test, the following parameters are considered:

| Parameter                          | Symbol     | Value                       |
| ---------------------------------- | ---------- | --------------------------- |
| Diffusion coefficient                | $M$   | $1.e-8$                       |
| Domain size              | $L$   | $1.e-3$                     |


### __Numerical scheme__

- Time marching: Euler Implicit scheme, $t\in[0,0.5]$, $\delta t=0.1$
- Spatial discretization: built from MFEM + FE order 1 

 



## Running

### __Using the binary__
```shell
./Diffusion1Dtest1
```

### __Using ctest__

```shell
ctest -R Diffusion1Dtest1
```

### __In case of code coverage analysis__

```shell
make Diffusion1Dtest1_coverage
```


## Post-processing

TO BE DONE 

## Files & Dependencies

Source file : `main.cpp`

## References

None

## Intellectual Property

See [About page](../../../../../about.html) 