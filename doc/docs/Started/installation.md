

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


## Installing SLOTH On Supercomputers Without Internet Access

This guide provides detailed steps to install SLOTH on a supercomputer without internet access using provided scripts. It includes using Spack for package management and compiling dependencies required by SLOTH.

Installing SLOTH on a supercomputer without internet access involves preparing the environment, downloading necessary components, creating a local Spack mirror, and building SLOTH with all dependencies.

### Use And Adapt Scripts

The installation is performed in two main parts using the scripts provided below:

1. `sloth-topaze-part1.sh`: Prepares the environment and creates an archive of necessary components.
2. `sloth-topaze-part2.sh`: Sets up Spack and compiles SLOTH on the target supercomputer.

Make sure to adapt the environment variables in the scripts (e.g., `MY_LOG`, `DEST_DIR`) to your specific user settings.

On your local machine:

```
source sloth-topaze-part1.sh
```

On your distant machine (Topaze in our example)

```
sloth-topaze-part2.sh
```

### Part 1: Preparing the Environment (`sloth-topaze-part1.sh`)

This script is designed to be run on a local machine with internet access. It sets up the environment, clones necessary repositories, prepares Spack, and packages everything into an archive for transfer to the supercomputer.

#### Step-by-Step Breakdown

1. **Define Root and Working Directories:**
   ```bash
   export ROOT_DIR=$PWD
   mkdir -p sloth-topaze-dir && cd sloth-topaze-dir
   export WORK_DIR=$ROOT_DIR/sloth-topaze-dir
   export MY_LOG=your_login      # Replace with your Topaze login
   export DEST_DIR=/path/to/destination # Replace with your destination directory
   ```
   - `ROOT_DIR` is set to the current directory.
   - Creates a subdirectory `sloth-topaze-dir` where all operations will occur.
   - `WORK_DIR` is set to the path of `sloth-topaze-dir`.
   - `MY_LOG` and `DEST_DIR` are placeholders for your supercomputer login and destination directory. You need to replace these with your actual login and path on the supercomputer.

2. **Clone Spack Repository:**
   ```bash
   echo "Getting Spack ..."
   if [ ! -d "spack" ]; then
       git clone https://github.com/spack/spack.git
   fi
   export SPACK_ROOT=$PWD/spack
   rm -r ~/.spack
   source ${SPACK_ROOT}/share/spack/setup-env.sh
   ```
   - Clones Spack from GitHub.
   - Sets `SPACK_ROOT` to the path of the cloned Spack directory.
   - Removes any existing `.spack` configuration to ensure a clean setup.
   - Sources the Spack environment to set up paths and commands for use.

3. **Clone SLOTH Repository:**
   ```bash
   echo "Getting PLEIADES/SLOTH ..."
   if [ ! -d "sloth" ]; then
       git clone https://www-git-cad.intra.cea.fr/DEC/collaboratif/ci230846/DEV_PROJECT/sloth.git
   fi
   ```
   - Similar to Spack, this step clones the SLOTH repository if it doesn't already exist in the working directory.

4. **Create a Spack Bootstrap Mirror:**
   ```bash
   spack bootstrap mirror --binary-packages my_bootstrap
   ```
   - Creates a bootstrap mirror that includes binary packages of the basic build tools that Spack needs to work offline.

5. **Create a Specific Spack Mirror for Dependencies:**
   ```bash
   spack mirror create -d mirror-mfem -D gcc@11.2.0 mfem+mpi+debug+openmp+petsc+strumpack+suite-sparse+sundials+superlu-dist%gcc@11.2.0
   ```
   - Creates a mirror named `mirror-mfem` for all specified dependencies (`mfem`, `petsc`, etc.), ensuring that Spack can access these packages without internet access on the supercomputer.
   - You can add extra packages here.

6. **Package and Transfer Files:**
   ```bash
   cd $ROOT_DIR
   tar cvf archive.tar.gz sloth-topaze-dir/
   scp archive.tar.gz $MY_LOG@topaze.ccc.cea.fr:$DEST_DIR/
   ```
   - Archives the entire `sloth-topaze-dir` directory into `archive.tar.gz`.
   - Uses `scp` to securely copy this archive to the specified destination directory on the supercomputer. Replace `topaze.ccc.cea.fr` with the appropriate hostname if needed.

---

### Part 2: Setting Up and Building SLOTH (`sloth-topaze-part2.sh`)

This script is run on the supercomputer. It unpacks the archive, sets up the Spack environment, configures Spack to work offline, and builds SLOTH with all required dependencies.

#### Step-by-Step Breakdown

1. **Define Directories:**
   ```bash
   export DEST_DIR=$PWD
   export WORK_DIR=$DEST_DIR/sloth-topaze-dir
   ```
   - `DEST_DIR` is set to the current working directory (where the archive was transferred).
   - `WORK_DIR` points to the `sloth-topaze-dir` directory inside `DEST_DIR`.

2. **Clean Up and Extract the Archive:**
   ```bash
   rm -r ~/.spack
   cd $DEST_DIR
   tar xvf archive.tar.gz
   cd $WORK_DIR
   ```
   - Removes any existing Spack configuration (`~/.spack`) to ensure a fresh environment setup.
   - Extracts the archive (`archive.tar.gz`) containing all previously prepared files.

3. **Set Up Spack Environment:**
   ```bash
   source $WORK_DIR/spack/share/spack/setup-env.sh
   spack bootstrap reset -y
   spack bootstrap add --scope=site --trust local-binaries $PWD/my_bootstrap/metadata/binaries/
   spack bootstrap disable --scope=site github-actions-v0.5
   spack bootstrap disable --scope=site github-actions-v0.4
   spack bootstrap disable --scope=site spack-install
   spack bootstrap root $PWD/spack/bootstrap
   spack bootstrap now
   spack bootstrap status
   ```
   - Sources the Spack environment to set up the paths and commands.
   - Resets Spack’s bootstrap configuration and adds the local bootstrap mirror (`my_bootstrap`) created earlier, ensuring all dependencies are fetched locally.
   - Disables unnecessary bootstrap sources (`github-actions-v0.5`, etc.) to avoid any attempt to connect online.
   - Sets the root path for Spack’s bootstrap environment and checks the status.

4. **Set Compiler and Environment Variables:**
   ```bash
   export CC='gcc'
   export CXX='g++'
   export FC='mpifort'
   export OMPI_CC='gcc'
   export OMPI_CXX='g++'
   export OMPI_FC='gfortran'
   ```
   - Specifies compilers for C, C++, and Fortran, ensuring the correct toolchain is used during the build.
   - Sets OpenMPI environment variables to link the compilers correctly.

5. **Load Required Modules and Add Spack Mirror:**
   ```bash
   module load gnu/11.2.0 mpi cmake/3.29.6 
   spack mirror add SLOTH $WORK_DIR/mirror-mfem/
   spack compiler find
   spack external find openmpi
   spack external find cmake
   spack external find openssh
   ```
   - Loads necessary modules (`gnu`, `mpi`, `cmake`) to provide the required tools and compilers.
   - Adds the previously created Spack mirror (`mirror-mfem`) so that dependencies are fetched from the local mirror instead of the internet.
   - Detects and registers available compilers and external software (e.g., `openmpi`, `cmake`, `openssh`) within Spack.

6. **Install Dependencies and Build SLOTH:**
   ```bash
   spack install gcc@11.2.0 mfem+mpi+debug+openmp+petsc+strumpack+suite-sparse+sundials+superlu-dist%gcc@11.2.0
   cd $WORK_DIR/sloth
   mkdir build && cd build
   spack load mfem
   spack load metis
   export HYPRE_DIR=`spack location -i hypre`
   export MPI_DIR=`spack location -i mpi`
   export METIS_DIR=`spack location -i metis`

   cmake .. -DMFEM_USE_PETSC=ON -DPETSC_DIR=${PETSC_DIR} -DPETSC_ARCH="" -DPETSC_INCLUDES=${PETSC_DIR}/include -DPETSC_LIBRARIES=${PETSC_DIR}/lib -DPETSC_EXECUTABLE_RUNS=${PETSC_DIR}/bin
   make -j 10
   ctest
   ```
   - Installs GCC and other specified dependencies from the local mirror without accessing the internet.
   - Sets up the build directory within the SLOTH repository (`build`).
   - Loads required dependencies (`mfem`, `metis`) to ensure they are available for the build process.
   - Sets environment variables to locate specific dependency installations.
   - Configures SLOTH with `cmake`, pointing to relevant dependencies (`PETSC`, etc.), and builds the software using `make`.
   - Runs tests with `ctest` to verify the build.

---

### Run Your Simulation On Topaze


Script example of a simulation running on milan partition over 8192 mpi processes with a duration limit of about 24 hours:

```
#!/bin/bash
#MSUB -r sloth_big_run
#MSUB -n 8192
#MSUB -c 1
#MSUB -T 86000
#MSUB -m scratch
#MSUB -o sloth_big_run_%I.o
#MSUB -e sloth_big_run_%I.e
#MSUB -q milan

set -x
export OMP_NUM_THREADS=1
ccc_mprun ./test3D

```

