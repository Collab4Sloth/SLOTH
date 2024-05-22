spack load mfem
export HYPRE_DIR=`spack location -i hypre`
export MPI_DIR=`spack location -i mpi`

cmake .. -DMFEM_USE_PETSC=ON -DPETSC_DIR=${PETSC_DIR} -DPETSC_ARCH="" -DPETSC_INCLUDES=${PETSC_DIR}/include -DPETSC_LIBRARIES=${PETSC_DIR}/lib -DPETSC_EXECUTABLE_RUNS=${PETSC_DIR}/bin -DCMAKE_BUILD_TYPE=$1