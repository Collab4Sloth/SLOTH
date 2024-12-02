# 1D simulation (temporary)


## Statement of the problem

 
### __Governing equation__
Let us consider the following set of governing equations:
$$
\frac{\partial c}{\partial t}=[\nabla \cdot{} (M \nabla c - M \nabla c)]\text{ in }\Omega 
$$



### __Boundary conditions__

Neumann boundary conditions are prescribed on $\Gamma$:
$$
{\bf{n}} \cdot{} \lambda \nabla c=0 \text{ on }\Gamma
$$


### __Parameters__
For this test, the following parameters are considered:

| Parameter                          | Symbol     | Value                       |
| ---------------------------------- | ---------- | --------------------------- |
| Primary coefficient                | $M$        | $1.e-8$                       |
| Domain size                        | $L$        | $1.e-3$                     |


### __Numerical scheme__

- Time marching: Euler Implicit scheme, $t\in[0,0.5]$, $\delta t=0.1$
- Spatial discretization: built from MFEM + order 1 FE 


 



## Running

### __Using the binary__
```shell
./Diffusion2Dtest1
```

### __Using ctest__

```shell
ctest -R Diffusion1Dtest3
```

### __In case of code coverage analysis__

```shell
make Diffusion1Dtest3_coverage
```

## Files & Dependencies

Source file : `main.cpp`

## References

None

## Intellectual Property

See [About page](../../../../../about.html) 