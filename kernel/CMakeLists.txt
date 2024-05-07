# Create the library

add_library(Sloth INTERFACE)

file(GLOB_RECURSE SRC_FILES ${CMAKE_SOURCE_DIR}/kernel/*.hpp)
target_sources(Sloth INTERFACE ${SRC_FILES})

target_compile_features(Sloth INTERFACE cxx_std_17)
set_target_properties(Sloth PROPERTIES LINKER_LANGUAGE CXX)

# And the most important part set the include directories
target_include_directories(Sloth INTERFACE ${CMAKE_CURRENT_LIST_DIR})

target_include_directories(Sloth INTERFACE ${MFEM_DIR}/include)
target_link_directories(Sloth INTERFACE ${MFEM_DIR}/lib $ENV{HYPRE_DIR}/lib $ENV{MPI_DIR}/lib $ENV{PETSC_DIR}/lib)
target_link_libraries(Sloth INTERFACE mfem mpi z HYPRE SuiteSparse::suitesparse metis petsc $ENV{FILESYSTEM_VAR})
