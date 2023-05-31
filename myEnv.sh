#!/bin/bash
# SPACK
export SPACK_ROOT=$HOME/home-local/MyGitProjects/spack
source ${SPACK_ROOT}/share/spack/setup-env.sh
# spack load hypre metis mgis@master cmake@3.21.4 mumps tfel@4.0.0 suite-sparse gcc@11.2.0 
spack load mfem  mpi  petsc strumpack suite-sparse sundials superlu-dist
# PUGIXML 
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

