#!/bin/sh
export DEST_DIR=$PWD
export WORK_DIR=$DEST_DIR/sloth-topaze-dir
rm -r ~/.spack
cd $DEST_DIR
tar xvf archive.tar.gz
cd $WORK_DIR
source $WORK_DIR/spack/share/spack/setup-env.sh
spack bootstrap reset -y
spack bootstrap add --scope=site --trust local-binaries $PWD/my_bootstrap/metadata/binaries/
spack bootstrap disable --scope=site github-actions-v0.5
spack bootstrap disable --scope=site github-actions-v0.4
spack bootstrap disable --scope=site spack-install
spack bootstrap root $PWD/spack/bootstrap
spack bootstrap now
spack bootstrap status
export CC='gcc'
export CXX='g++'
export FC='mpifort'
export OMPI_CC='gcc'
export OMPI_CXX='g++'
export OMPI_FC='gfortran'

module load gnu/11.2.0 mpi cmake/3.29.6 
spack mirror add SLOTH $WORK_DIR/mirror-mfem/
spack compiler find find 
spack external find openmpi
spack external find cmake
spack external find openssh
spack external find gmake
spack install mfem+mpi+openmp+petsc+strumpack+suite-sparse+sundials+superlu-dist
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

