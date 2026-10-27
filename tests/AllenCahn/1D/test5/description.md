# 1D Allen-Cahn simulation with adaptive time-stepping: rejected time-steps

Same physical case as `AllenCahn1Dtest3`, $t\in[0,100]$.

The Newton solver is limited to 3 iterations: the first time-step ($\delta t = 50$) fails. It is rejected, the variables (unknowns and history) are restored, and the time-step is retried with a smaller time-step. The time-step then grows again when Newton converges in at most 2 iterations.

`newton_iterations_high` is greater than `iter_max`: the time-step is reduced by rejections only. The comparison with the reference `time_specialized.csv` checks the restoration of the variables after a rejection.

| Parameter                  | Value |
| -------------------------- | ----- |
| `time_step` (initial)      | 50    |
| `time_step_min`            | 1.e-6 |
| `time_step_cutback_factor` | 0.5   |
| `time_step_growth_factor`  | 1.5   |
| `newton_iterations_low`    | 2     |
| `newton_iterations_high`   | 4     |
| `max_time_step_rejections` | 30    |
| Newton `iter_max`          | 3     |

## Running

```shell
ctest -R AllenCahn1Dtest5
```

## Files & Dependencies

Source file : `main.cpp`

## Intellectual Property

See [About page](../../../../../about.html)
