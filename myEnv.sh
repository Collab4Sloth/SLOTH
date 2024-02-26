#!/bin/bash
# SPACK
export SPACK_ROOT=$HOME/home-local/MyGitProjects/spack
source ${SPACK_ROOT}/share/spack/setup-env.sh
spack load mfem mpi petsc strumpack suite-sparse sundials superlu-dist
export SPACK=${SPACK_ROOT}/opt/spack/linux-debian10-skylake/gcc-8.3.0/
# DEPENDENCES
export OPENMPI=`spack location -i openmpi`
export HYPRE=`spack location -i hypre`
export SUPERLUDIST=`spack location -i superlu-dist`
export PARMETIS=`spack location -i parmetis`
export METIS=`spack location -i metis`
export SUNDIALS=`spack location -i sundials`
export SUITE=`spack location -i suite-sparse`
export STRUM=`spack location -i strumpack`
export NETLIB=`spack location -i netlib-scalapack`
export BUTTERFLY=`spack location -i butterflypack`
export ZFP=`spack location -i zfp`
export PETSC=`spack location -i petsc`
export ZLIB=`spack location -i zlib`

# MFEM

export MFEM_INSTALL_DIR=`spack location -i mfem`
# export PUGIXML_INSTALL_DIR=`spack location -i pugixml`
export MFEM_INCLUDE=$MFEM_INSTALL_DIR/include
export MFEM_INCLUDE_MFEM=$MFEM_INCLUDE/mfem
export MFEM_INCLUDE_MFEM_FEM=$MFEM_INCLUDE_MFEM/fem
export MFEM_INCLUDE_MFEM_GENERAL=$MFEM_INCLUDE_MFEM/general
export MFEM_INCLUDE_MFEM_LINALG=$MFEM_INCLUDE_MFEM/linalg
export MFEM_INCLUDE_MFEM_MESH=$MFEM_INCLUDE_MFEM/mesh
# more $MFEM_INSTALL_DIR/share/mfem/config.mk
export SRC_DIR=$PWD

