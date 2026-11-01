# 1D Allen-Cahn simulation with adaptive time-stepping


## Statement of the problem
This example code solves one Allen-Cahn equation in a 1D domain $\Omega=[0,R]$ with $R=10^{-3}$ (inline mesh), as in `AllenCahn1Dtest1`, with an adaptive time-step.

The initial time-step covers the whole time interval and the number of Newton iterations is limited: the first attempts fail. Each failed time-step is rejected (the variables are restored) and retried with a smaller time-step. The time-step then grows again when Newton converges in few iterations.

### __Governing equation__

$$
\frac{\partial \phi}{\partial t}=M_\phi[\nabla \cdot{} \lambda \nabla \phi -\omega W'(\phi)]\text{ in }\Omega
$$

$$
W(\phi)=\phi^2(1-\phi)^2
$$

### __Boundary conditions__

Neumann boundary conditions are prescribed on $\Gamma_{left}$ (r=0) and $\Gamma_{right}$ (r=R).

### __Parameters__

| Parameter                          | Symbol     | Value                       |
| ---------------------------------- | ---------- | --------------------------- |
| mobility coefficient               | $M_\phi$   | $1.e-3$                     |
| energy gradient coefficient        | $\lambda$  | $\frac{3}{2}\sigma\epsilon$ |
| surface tension                    | $\sigma$   | $0.06$                      |
| interface thickness                | $\epsilon$ | $5.e-4$                     |
| depth of the double-well potential | $\omega$   | $12{\sigma}/{\epsilon}$     |


### __Numerical scheme__

- Time marching: Euler Implicit scheme, $t\in[0,100]$
- Adaptive time-stepping (`adaptive_time_step = true`):

| Parameter                  | Value |
| -------------------------- | ----- |
| `time_step` (initial)      | 50    |
| `time_step_min`            | 1.e-6 |
| `time_step_cutback_factor` | 0.5   |
| `time_step_growth_factor`  | 1.5   |
| `newton_iterations_low`    | 2     |
| `newton_iterations_high`   | 3     |
| `max_time_step_rejections` | 30    |

- Newton solver: at most 5 iterations, relative tolerance $10^{-10}$
- Spatial discretization: uniform mesh with $N=30$ elements

`AllenCahn1Dtest4` solves the same case with `max_time_step_rejections = 0`: the rejection of the first time-step stops the computation (expected failure).


## Running

### __Using the binary__
```shell
mpirun -np 2 ./AllenCahn1Dtest3
```

### __Using ctest__

```shell
ctest -L ADAPTIVE_DT
```


## Files & Dependencies

Source file : `main.cpp`

## References

None

## Intellectual Property

See [About page](../../../../../about.html)
