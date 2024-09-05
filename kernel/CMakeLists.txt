# Create the library

add_library(Sloth INTERFACE)

file(GLOB_RECURSE SRC_FILES ${CMAKE_SOURCE_DIR}/kernel/sloth.hpp 
                            ${CMAKE_SOURCE_DIR}/kernel/BCs/*.hpp               
                            ${CMAKE_SOURCE_DIR}/kernel/Coefficients/*.hpp          
                            ${CMAKE_SOURCE_DIR}/kernel/Convergence/*.hpp    
                            ${CMAKE_SOURCE_DIR}/kernel/Couplings/*.hpp  
                            ${CMAKE_SOURCE_DIR}/kernel/Integrators/*.hpp         
                            ${CMAKE_SOURCE_DIR}/kernel/Operators/*.hpp         
                            ${CMAKE_SOURCE_DIR}/kernel/Parameters/*.hpp   
                            ${CMAKE_SOURCE_DIR}/kernel/PostProcessing/*.hpp      
                            ${CMAKE_SOURCE_DIR}/kernel/Problems/*.hpp        
                            ${CMAKE_SOURCE_DIR}/kernel/Profiling/*.hpp      
                            ${CMAKE_SOURCE_DIR}/kernel/Solvers/*.hpp          
                            ${CMAKE_SOURCE_DIR}/kernel/Spatial/*.hpp        
                            ${CMAKE_SOURCE_DIR}/kernel/Time/*.hpp               
                            ${CMAKE_SOURCE_DIR}/kernel/Utils/*.hpp             
                            ${CMAKE_SOURCE_DIR}/kernel/Variables/*.hpp  )

target_sources(Sloth INTERFACE ${SRC_FILES})


target_compile_features(Sloth INTERFACE cxx_std_17)
set_target_properties(Sloth PROPERTIES LINKER_LANGUAGE CXX)

# And the most important part set the include directories
target_include_directories(Sloth INTERFACE ${CMAKE_CURRENT_LIST_DIR})

target_include_directories(Sloth SYSTEM INTERFACE ${MFEM_DIR}/include)
target_link_directories(Sloth INTERFACE ${MFEM_DIR}/lib  $ENV{HYPRE_DIR}/lib $ENV{MPI_DIR}/lib $ENV{PETSC_DIR}/lib $ENV{METIS_DIR}/lib)

target_link_libraries(Sloth INTERFACE mfem mpi z HYPRE SuiteSparse::suitesparse metis $ENV{FILESYSTEM_VAR})

# PETSC ??
if(MFEM_USE_PETSC)
    target_link_libraries(Sloth INTERFACE petsc )
endif(MFEM_USE_PETSC)
