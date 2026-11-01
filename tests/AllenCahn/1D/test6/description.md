# 1D Allen-Cahn simulation with adaptive time-stepping: BDF2 scheme

Same case as `AllenCahn1Dtest5` ($t\in[0,100]$, Newton limited to 3 iterations) with the variable-step BDF2 scheme.

- The first time-steps fail (Newton not converged): they are rejected. The variables, their history (which gives $u^{n-1}$ to BDF2), their number of updates and the previous time-step $\delta t_{n-1}$ of BDF2 are restored, so that the time-step is retried with the ratio $\omega = \delta t_n/\delta t_{n-1}$ of the last accepted time-step.
- The growth factor (2.5) is larger than the ratio allowed between two successive time-steps with BDF2 (`time_step_max_ratio`, 2 by default, lower than $1+\sqrt{2}$ for zero-stability): each increase of the time-step is limited to a factor 2.

| Parameter                  | Value |
| -------------------------- | ----- |
| Time scheme                | BDF2  |
| `time_step` (initial)      | 50    |
| `time_step_min`            | 1.e-6 |
| `time_step_cutback_factor` | 0.5   |
| `time_step_growth_factor`  | 2.5   |
| `newton_iterations_low`    | 2     |
| `newton_iterations_high`   | 4     |
| `max_time_step_rejections` | 30    |
| Newton `iter_max`          | 3     |

## Running

```shell
ctest -R AllenCahn1Dtest6
```

## Files & Dependencies

Source file : `main.cpp`

## Intellectual Property

See [About page](../../../../../about.html)
