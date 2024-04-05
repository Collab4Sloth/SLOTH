cmake_minimum_required(VERSION 3.8)
project(sloth)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/INSTALL)

# Check GCC version
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)

  if(GCC_VERSION VERSION_LESS 9)
    find_package(Filesystem REQUIRED)
    set(ENV{FILESYSTEM_VAR} std::filesystem)
  else()
    set(ENV{FILESYSTEM_VAR})
  endif()
endif()

# (must be placed *before* any add_subdirectory, cmake bug ?)
include(CTest)
enable_testing()

# MFEM
find_package(MFEM REQUIRED)

if(MFEM_USE_MPI)
  message(STATUS "MFEM: using MPI, so SLOTH will use MPI (with HYPRE, METIS)")

  # The following include directories are automatically filled within FindMFEM.cmake
  include_directories(SYSTEM ${HYPRE_INCLUDE_DIRS})
  include_directories(SYSTEM ${METIS_INCLUDE_DIRS})
endif(MFEM_USE_MPI)

# Support for OpenMP
if(MFEM_USE_OPENMP OR MFEM_USE_LEGACY_OPENMP)
  message(STATUS "MFEM: using OPENMP, so SLOTH will use OPENMP")
  find_package(OpenMP REQUIRED)

  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif(MFEM_USE_OPENMP OR MFEM_USE_LEGACY_OPENMP)

if(MFEM_USE_MUMPS)
  find_package(OpenMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(MFEM_USE_MUMPS)

# Petsc
if(MFEM_USE_PETSC)
  message(STATUS "MFEM: using PETSc, so SLOTH will use PETSc")
  find_package(PETSc REQUIRED)
endif(MFEM_USE_PETSC)

# documentation
if(NOT TARGET doc)
  add_custom_target(doc)
  add_subdirectory(doc)
endif(NOT TARGET doc)

if(MFEM_USE_SUITESPARSE)
  find_package(SuiteSparse REQUIRED)
endif(MFEM_USE_SUITESPARSE)

find_package(METIS REQUIRED)

# HYPRE version
include(cmake/modules/hypre.cmake)
set(MFEM_HYPRE_VERSION ${HYPRE_VERSION})
add_compile_definitions(MFEM_HYPRE_VERSION=${MFEM_HYPRE_VERSION})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(Tests)
