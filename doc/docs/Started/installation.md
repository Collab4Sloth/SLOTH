

The installation of SLOTH consists of installing MFEM first and then, to compile SLOTH 

### Installation of MFEM on Linux using spack

A straightforward way to install MFEM is to use [spack](https://spack.readthedocs.io/en/latest/getting_started.html)

- First, clone spack and install it into `$SPACK` directory (see [spack](https://spack.readthedocs.io/en/latest/getting_started.html))
- Second, run the following commands to install mfem with right additional packages

```shell
$SPACK/share/spack/setup-env.sh

spack install mfem+mpi+debug+openmp+petsc+strumpack+suite-sparse+sundials+superlu-dist
```

____

### Installation of MFEM on Mac Os using homebrew 

MFEM can be install using homebrew. 

```shell
brew install mfem
```

by default, this installation depends on hypre, metis, openblas, suite-sparse.

It is possible rebuild  mfem with additional dependencies. 

- Get the .rb file : run `brew edit mfem` to open the default rb file or get it from [Github](https://github.com/Homebrew/homebrew-core/blob/5ecde7427aa47ac931795c78669f0a4da53a12ed/Formula/m/mfem.rb)
- Add your dependencies with `depends_on` directive. Here, let us consider the `petsc` dependency:

```shell
  depends_on "cmake" => :build
  depends_on "hypre"       
  depends_on "metis"       
  depends_on "openblas"
  depends_on "suite-sparse"
  depends_on "petsc"
```

- save the file in the directory and run the following command:
  
```shell
  brew install --formula mfem.rb
````

Installation with petsc can be checked by editing once again the mfem.rb file. petsc must be mentioned as default dependency. 

Each dependency can be installed easily using homebrew. 

## Compiling SLOTH using mfem version installed with homebrew

```shell
mkdir build ; cd build
```


- Set several environment variables 

```shell
export MFEM_DIR=$(echo `brew --prefix mfem`)

export MPI_DIR=$(echo `brew --prefix open-mpi`)

export HYPRE_DIR=$(echo `brew --prefix hypre`)

export METIS_DIR=$(echo `brew --prefix metis`)

```

____

### Compiling SLOTH


- First, create a dedicated directory to build SLOTH
```shell
mkdir build && cd build
```

- Second, load the SLOTH environment file 

```shell
source ../envSloth.sh [OPTIONS] 
```

where [OPTIONS] are:
```shell
    --release to build with Release compiler options 
        
    --debug to build with Debug compiler options 
        
    --coverage to build with Coverage compiler options 
        
    --clean to remove previous build if it exists 
```

By default, SLOTH is built with release compiler options.


- Finally, compile 

```shell
make -j N 
```
with N the number of jobs.

## Static code analysis

- To do a static code analysis with CPPLINT, run 

```shell
make lint
```

## Building documentation

- To build documentation, run 

```shell
make doc
```


## Building/running tests

- All tests are located in the Tests subdirectory. 
- All tests are built during compilation (see `make` command). 
- All tests can be executed by running
```shell
make test 
```

or 

```shell
ctest
```

- This stage can be mimicked without execution by using the `-N` option:

```shell
ctest -N
```

- Each test can be run independently by using the `-R` option:
- 
```shell
ctest -R monTest
```


The documentation is built with [MKdocs](https://www.mkdocs.org/) and, for code documentation, with Doxygen. 

[MKdocs](https://www.mkdocs.org/) is used to generate a global documentation under HTML format. A PDF file is also generated. 
To do it, a set of packages are needed. Please refer to the following list (may be not complete):
```shell

pip3 install mkdocs                           # for the use of MkDocs
pip3 install mkdocs-material                  # for the graphical environment Material
pip3 install mkdocs-material-extensions       # for additional functionalities of Material
pip3 install mkdocs-plugin-offline            # for the plugin offline
pip3 install mkdocs-plugin-search             # for the plugin search
pip3 install mkdocs-with-pdf                  # for the plugin with-pdf

pip3 install markdown                         # for the extensions attr_list, md_in_html, def_list
pip3 install pymdown-extensions               # for the extensions highlight, inlinehilite, snippets, superfences, tasklist, arithmatex



pip3 install  mkdocs-git-revision-date-plugin  # to include Git revision....
pip3 install  mkdocs-material[imaging]         # to use additional graphical features of  Material.
pip3 install  mkdocs-pdf-export-plugin         # To export documentation under PDF format
pip3 install  python-markdown-math             # To take into account mathematical formula in markdown
pip3 install  WeasyPrint==52.5                 # Used by the plugin with-pdf.


```





