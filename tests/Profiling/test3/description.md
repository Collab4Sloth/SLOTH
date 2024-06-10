# Profiling test 3

## Statement of the problem
This example code measures the execution time of Profiling functions through a test code. It begins by initiating global profiling, within which the profiling for the test code is activated. After finishing of the global profiling, The difference in time between the two profiles is calculated, representing the execution time of Profiling functions.

## Running

### __Using the binary__
```shell
./Profilingtest3
```

### __Using ctest__

```shell
ctest -R Profilingtest3
```

### __In case of code coverage analysis__

```shell
make Profilingtest3_coverage
```