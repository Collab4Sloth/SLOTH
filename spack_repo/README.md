# Spack repository for SLOTH

This directory is a [Spack](https://spack.readthedocs.io) package repository
containing the `sloth` package. Spack installs MFEM, its dependencies and SLOTH.

## Install SLOTH

```bash
source $SPACK/share/spack/setup-env.sh
spack repo add /path/to/SLOTH/spack_repo
spack install sloth
```

Versions and variants:

```bash
spack info sloth
```

| Variant | Default | Effect |
|---|---|---|
| `+shared` | off | shared library (requires an MFEM built with `+shared`) |
| `+mumps` | off | MUMPS direct solver, through MFEM |
| `+voro` | off | Voro++ (Voronoi tessellations) |
| `+bigint` | off | 64-bit global indices (HYPRE `+mixedint`), for more than ~2 billion unknowns |
| `+native` | off | `-march=native` and link-time optimisation (not portable to other CPUs) |

Examples:

```bash
spack install sloth@master            # development version
spack install sloth+mumps
spack install sloth %gcc@14           # another compiler
```

An MFEM already installed with Spack is reused when it satisfies the requirements
of the package (`mfem@4.9: +mpi`).

## Use the installed SLOTH

```bash
spack load sloth
cp -r $(spack location -i sloth)/share/sloth/examples/minimal .
cd minimal
cmake -B build && cmake --build build
cd build && ctest --output-on-failure
```

`spack load sloth` makes SLOTH visible to CMake (`find_package(Sloth)`) and provides
Python with SymPy for the coefficient generator. See the example's `CMakeLists.txt`
and the documentation page "Using an installed SLOTH in your own project".
