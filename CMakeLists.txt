cmake_minimum_required(VERSION 3.8)
project(sloth)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# (must be placed *before* any add_subdirectory, cmake bug ?)
include(CTest)
enable_testing()

# MFEM
find_package(MFEM REQUIRED)

if (MFEM_USE_MPI)
  message(STATUS "MFEM: using MPI, so SLOTH will use MPI")
  find_package(MPI REQUIRED)
  # The following include directories are automatically filled within FindMFEM.cmake
  include_directories(SYSTEM ${HYPRE_INCLUDE_DIRS})
  include_directories(SYSTEM ${METIS_INCLUDE_DIRS})
endif(MFEM_USE_MPI)

if (MFEM_USE_MUMPS)
	find_package(OpenMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif (MFEM_USE_MUMPS)

# Support for OpenMP
if(MFEM_USE_OPENMP OR MFEM_USE_LEGACY_OPENMP) 
  find_package(OpenMP REQUIRED)
endif(MFEM_USE_OPENMP OR MFEM_USE_LEGACY_OPENMP)
 
# documentation
if (NOT TARGET doc)
  add_custom_target(doc)
endif(NOT TARGET doc)

find_package(SuiteSparse REQUIRED)
find_package(METIS REQUIRED)

# HYPRE version
include(cmake/modules/hypre.cmake)  
set(MFEM_HYPRE_VERSION ${HYPRE_VERSION})
add_compile_definitions(MFEM_HYPRE_VERSION=${MFEM_HYPRE_VERSION})

message($ENV{HYPRE_DIR})
message(${MFEM_DIR})

message("INCLUDE_DIR " ${CMAKE_CURRENT_SOURCE_DIR})
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
add_subdirectory(Tests)

