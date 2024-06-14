

export PETSC_DIR=`echo $(brew --prefix petsc)` 
export HYPRE_DIR=`echo $(brew --prefix hypre)`
export MPI_DIR=`echo $(brew --prefix open-mpi)`
export METIS_DIR=`echo $(brew --prefix metis)`
export MFEM_DIR=$(echo `brew --prefix mfem`)

cmake .. -DMFEM_USE_PETSC=ON -DPETSC_DIR=${PETSC_DIR} -DPETSC_ARCH=""  -DCMAKE_BUILD_TYPE=$1