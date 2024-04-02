# How To install Sloth with cmake and spack

```
mkdir build && cd build
spack load mfem
spack load hypre
export HYPRE_DIR=`spack location -i hypre`
cmake ..
make

./Tests/AllenCahn/1D/test1/AllenCahn1Dtest1 
./Tests/AllenCahn/3D/test1/AllenCahn3Dtest1 
```
