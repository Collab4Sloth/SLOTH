# Create the library

add_library(Sloth INTERFACE)

file(GLOB_RECURSE SRC_FILES ${CMAKE_SOURCE_DIR}/kernel/sloth.hpp 
                            ${CMAKE_SOURCE_DIR}/kernel/BCs/*.hpp               
                            ${CMAKE_SOURCE_DIR}/kernel/Coefficients/*.hpp          
                            ${CMAKE_SOURCE_DIR}/kernel/Convergence/*.hpp    
                            ${CMAKE_SOURCE_DIR}/kernel/Couplings/*.hpp  
                            ${CMAKE_SOURCE_DIR}/kernel/Integrators/*.hpp       
                            ${CMAKE_SOURCE_DIR}/kernel/Options/*.hpp           
                            ${CMAKE_SOURCE_DIR}/kernel/Operators/*.hpp         
                            ${CMAKE_SOURCE_DIR}/kernel/Parameters/*.hpp   
                            ${CMAKE_SOURCE_DIR}/kernel/PostProcessing/*.hpp      
                            ${CMAKE_SOURCE_DIR}/kernel/Problems/*.hpp        
                            ${CMAKE_SOURCE_DIR}/kernel/Profiling/*.hpp      
                            ${CMAKE_SOURCE_DIR}/kernel/Solvers/*.hpp          
                            ${CMAKE_SOURCE_DIR}/kernel/Spatial/*.hpp        
                            ${CMAKE_SOURCE_DIR}/kernel/Time/*.hpp               
                            ${CMAKE_SOURCE_DIR}/kernel/Utils/*.hpp             
                            ${CMAKE_SOURCE_DIR}/kernel/Variables/*.hpp         
                            ${CMAKE_SOURCE_DIR}/tests/tests.hpp  )

target_sources(Sloth INTERFACE ${SRC_FILES})


# target_compile_features(Sloth INTERFACE cxx_std_20)
set_target_properties(Sloth PROPERTIES LINKER_LANGUAGE CXX)

# And the most important part set the include directories
target_include_directories(Sloth INTERFACE ${CMAKE_CURRENT_LIST_DIR})

target_include_directories(Sloth SYSTEM INTERFACE ${MFEM_DIR}/include)
target_link_directories(Sloth INTERFACE ${MFEM_DIR}/lib  $ENV{HYPRE_DIR}/lib $ENV{MPI_DIR}/lib $ENV{PETSC_DIR}/lib $ENV{METIS_DIR}/lib )

target_link_libraries(Sloth INTERFACE ${MFEM_LIBRARY} ${MPI_CXX_LIBRARIES} z ${HYPRE_LIBRARY} ${SuiteSparse_LIBRARIES} ${METIS_LIBRARY} $ENV{FILESYSTEM_VAR})


# PETSC ??
if(MFEM_USE_PETSC)
    target_link_libraries(Sloth INTERFACE petsc )
endif(MFEM_USE_PETSC)


# EXTERNAL??
option(SLOTH_USE_EXTERNAL "Enable the use of an external package that requires linking with external library" OFF)
option(SLOTH_USE_AUTO_EXTERNAL "Enable the use of an autonomous external package" OFF)
if(SLOTH_USE_EXTERNAL OR SLOTH_USE_AUTO_EXTERNAL)
    file(GLOB_RECURSE EXT_SRC_FILES $ENV{EXT_SRC}/*.hpp)
    add_library(SlothExternal INTERFACE ${EXT_SRC_FILES})
    if(SLOTH_USE_EXTERNAL)
        target_link_directories(Sloth INTERFACE $ENV{EXT_LIBDIR} )
        target_link_libraries(Sloth INTERFACE SlothExternal $ENV{EXT_LIBNAME})
    endif(SLOTH_USE_EXTERNAL)
endif()


# LibTorch API ??
if(SLOTH_USE_LIBTORCH)
    find_package(Torch REQUIRED)
    MESSAGE(STATUS "Sloth use libTorch API ")
    target_include_directories(Sloth INTERFACE ${TORCH_INCLUDE_DIRS} )
    target_link_libraries(Sloth INTERFACE ${TORCH_LIBRARIES})
endif(SLOTH_USE_LIBTORCH)

