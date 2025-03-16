Compilation:
############

Get MFEM:

```
spack install mfem+mpi+suite-sparse+sundials+superlu-dist+shared
spack load mfem
export HYPRE_DIR=`spack location -i hypre`
```


Get Onika and compile it:

```
git clone https://github.com/Collab4exaNBody/onika.git
mkdir install-onika
export export onika_DIR=$PWD/install-onika
cd onika
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$onika_DIR -DCMAKE_BUILD_TYPE=Release
make install -j 10
```


Compile the SLOTH Proto:

```
https://github.com/Collab4Sloth/SLOTH.git
export ONIKA_CONFIG_PATH=$PWD/SLOTH/prototype/data/config/
export ONIKA_PLUGIN_PATH=$PWD/SLOTH/build/prototype/plugins/
mkdir SLOTH/build && cd SLOTH/build
cmake .. -DCMAKE_BUILD_TYPE=Release
cd prototype
make -j 4
```

Comments and remarks
####################

This is a list of small comments to tend the current SLOTH design to the SLOTH Proto.

- Need to decrease the level of template, especially for coupling and problem classes
- Add move operator, ie operator=(T&& value) for `*my_slot = std::move(data_to_keep_alive);`
- Split the `solve` function into several operators. The time discretization should be handled by Onika with a compute loop
- It's no longer required to define big classes that encompass several small classes (with pointers or references), it will be better to use Onika slots to get all you need. Note that onika can eliminate unneeded initialization.
- Need to add a definition of a default template for operators, I mean that every operator that needs to define dim or FEElement could be templatized and the user only needs to define an operator: `type_h1_dim_3` at the beginning of its input file.

