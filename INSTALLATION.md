# How To install Sloth with cmake and spack

The installation of SLOTH consists of installing MFEM first and then, to compile SLOTH 

### Installing MFEM

A straightforward way to install MFEM is to use [spack](https://spack.readthedocs.io/en/latest/getting_started.html)

- First, clone spack and install it into `$SPACK` directory (see [spack](https://spack.readthedocs.io/en/latest/getting_started.html))
- Second, run the following commands to install mfem with right additional packages

```shell
$SPACK/share/spack/setup-env.sh

spack install mfem+mpi+debug+openmp+petsc+strumpack+suite-sparse+sundials+superlu-dist
```

### Compiling SLOTH

- First, create a dedicated directory to build SLOTH
```shell
mkdir build && cd build
```

- Second, load mfem as SLOTH's prerequisite using `spack`
```shell
spack load mfem
```

- Third, export `HYPRE` and `MPI`  location into environment variables used during compilation process

```shell
export HYPRE_DIR=`spack location -i hypre`
export MPI_DIR=`spack location -i mpi`
```

- Fourth, run `cmake` with `PETSc` directives

```shell
cmake .. -DMFEM_USE_PETSC=ON -DPETSC_DIR=${PETSC_DIR} -DPETSC_ARCH="" -DPETSC_INCLUDES=${PETSC_DIR}/include -DPETSC_LIBRARIES=${PETSC_DIR}/lib -DPETSC_EXECUTABLE_RUNS=${PETSC_DIR}/bin
```
- Finally, compile (i.e. build all examples)

```shell
make