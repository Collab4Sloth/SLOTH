# 1D semi-infinite domain simulation


## Statement of the problem
This example code is based on the analytical solution for a semi-infinite domain.


 
### __Governing equation__
Let us consider the following set of governing equations:
$$
\frac{\partial c}{\partial t}=[\nabla \cdot M\nabla c]\text{ in }\Omega 
$$



### __Boundary conditions__

Neumann conditions are prescribed on boundaries:
$$
{\bf{n}} \cdot{} M \nabla c=0 \text{ on }\Gamma_{right} 
$$
$$
{\bf{n}} \cdot{} M \nabla c=0 \text{ on }\Gamma_{left}
$$

### __Initial conditions__

The composition is assumed to be constant at $t=0$: 
$$
c(x,t=0) = \frac{1}{2}\left[1+\tanh\left(\dfrac{x-L/2}{\epsilon}\right)\right] 
$$

### __Analytical solution__

The analytical solution can be written:
$$
c(x,t) = \frac{1}{2}\left[1+\text{erf}\left(\dfrac{x-L/2}{\sqrt{4Mt}}\right)\right] 
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
./Diffusion1Dtest2
```

### __Using ctest__

```shell
ctest -R Diffusion1Dtest2
```

### __In case of code coverage analysis__

```shell
make Diffusion1Dtest2_coverage
```


## Post-processing

TO BE DONE 

## Files & Dependencies

Source file : `main.cpp`

## References

None

## Intellectual Property

See [About page](../../../../../about.html) 