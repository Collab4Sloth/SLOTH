
set(CURRENT_EXE "AllenCahn2DPeriodictest2")
add_executable( ${CURRENT_EXE} main.cpp )
target_include_directories( ${CURRENT_EXE} PUBLIC ${MFEM_DIR}/include)
target_link_directories(${CURRENT_EXE} PUBLIC ${MFEM_DIR}/lib $ENV{HYPRE_DIR}/lib $ENV{MPI_DIR}/lib $ENV{PETSC_DIR}/lib)
target_link_libraries(${CURRENT_EXE} mfem mpi z HYPRE SuiteSparse::suitesparse metis petsc)
