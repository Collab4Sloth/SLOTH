SLOTH 
=============

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
```

## Building documentation

- To build documentation, run 

```shell
make doc
```

The documentation is built with [MKdocs](https://www.mkdocs.org/) and, for code documentation, with Doxygen. 

[MKdocs](https://www.mkdocs.org/) is used to generate a global documentation under HTML format. A PDF file is also generated. 
To do it, a set of packages are needed. Please refer to the following list (may be not complete):
```shell

pip3 install mkdocs                           # for the use of MkDocs
pip3 install mkdocs-material                  # for the graphical environment Material
pip3 install mkdocs-material-extensions       # for additional functionalities of Material
pip3 install mkdocs-with-pdf                  # for the plugin with-pdf

pip3 install markdown                         # for the extensions attr_list, md_in_html, def_list
pip3 install pymdown-extensions               # for the extensions highlight, inlinehilite, snippets, superfences, tasklist, arithmatex



pip3 install  mkdocs-git-revision-date-plugin  # to include Git revision....
pip3 install  mkdocs-material[imaging]         # to use additional graphical features of  Material.
pip3 install  mkdocs-pdf-export-plugin         # To export documentation under PDF format
pip3 install  python-markdown-math             # To take into account mathematical formula in markdown
pip3 install  WeasyPrint==52.5                 # Used by the plugin with-pdf.


```


<!-- ## Building documentation

To deploy the project, run 

```shell
make install 
``` -->

