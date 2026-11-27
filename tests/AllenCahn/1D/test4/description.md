# 1D Allen-Cahn simulation with adaptive time-stepping: no rejection allowed

Same case as `AllenCahn1Dtest3`, with `max_time_step_rejections = 0`.

The first time-step fails (Newton not converged): its rejection exceeds the maximum number of rejections and stops the computation. The test is registered as an expected failure: it checks that the first time-step of `AllenCahn1Dtest3` is actually rejected.

## Running

```shell
ctest -R AllenCahn1Dtest4
```

## Files & Dependencies

Source file : `main.cpp`

## Intellectual Property

See [About page](../../../../../about.html)
