# Profiling test 1

## Statement of the problem
This code example calculates the sum of numbers from 1 to n using MPI across multiple processes. Each process computes its local sum, which is then combined using MPI_Reduce to obtain the global sum. 

Additionally, this example includes profiling classes (Timers, Output, and UtilsforOutput) for performance analysis. Profiling can be enabled or disabled using the functions get_enableOutput() and get_disableOutput().

## Running

### __Using the binary__
```shell
mpirun -np <num_processes> ./Profilingtest1
```
Replace <num_processes> with the desired number of processes to run the code.

### __Using ctest__

```shell
ctest -R Profilingtest1
```

### __In case of code coverage analysis__

```shell
make Profilingtest1_coverage
```