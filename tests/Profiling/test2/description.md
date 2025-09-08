# Profiling test 2

## Statement of the problem
This MPI code prints a "Hello, World!" message from each process and calculates the elapsed time using MPI_Wtime and profiling classes in order to compare the diffenrece between the two methods.

Profiling can be enabled or disabled using the functions get_enableOutput() and get_disableOutput().

## Running

### __Using the binary__
```shell
mpirun -np <num_processes> ./Profilingtest2
```
Replace <num_processes> with the desired number of processes to run the code.

### __Using ctest__

```shell
ctest -R Profilingtest2
```

### __In case of code coverage analysis__

```shell
make Profilingtest2_coverage
```