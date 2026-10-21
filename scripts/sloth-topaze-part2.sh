#!/bin/bash
# Part 2 (on Topaze, no Internet access): build MFEM with Spack from the
# mirror of part 1, then build and test SLOTH.
export DEST_DIR=$PWD
export WORK_DIR=$DEST_DIR/sloth-topaze-dir

# Same MFEM spec in part 1 (mirror) and part 2 (installation)
MFEM_SPEC="mfem@4.9.0+mpi+mumps"

rm -rf ~/.spack
cd $DEST_DIR
tar xzvf archive.tar.gz
cd $WORK_DIR
source $WORK_DIR/spack/share/spack/setup-env.sh
# Package recipes copied in part 1 (no download from Topaze)
spack repo set --destination $WORK_DIR/spack-packages builtin

spack bootstrap reset -y
spack bootstrap add --scope=site --trust local-binaries $PWD/my_bootstrap/metadata/binaries/
# Check the names of the bootstrap sources with "spack bootstrap list"
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
spack compiler find
spack external find openmpi
spack external find cmake
spack external find openssh
spack external find gmake
spack install $MFEM_SPEC

# SLOTH: envSloth.sh finds MFEM and HYPRE with "spack location -i"
# (MFEM's third-party libraries are read from its config.mk).
# --portable: no -march=native, the login node may differ from compute nodes.
cd $WORK_DIR/sloth
mkdir -p build && cd build
bash ../envSloth.sh --portable --np=10
make -j 10 tests   # tests are not built by "make"
ctest
