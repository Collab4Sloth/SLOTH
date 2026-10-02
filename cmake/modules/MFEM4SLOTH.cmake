# MFEM
find_package(MFEM REQUIRED)

set(MFEM_INC "${MFEM_INSTALL_DIR}/include/mfem")

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

# MUMPS
if(MFEM_USE_MUMPS)
  message(STATUS "MFEM: using MUMPS, so SLOTH will use MUMPS")

  # MUMPS (spack)
  find_package(OpenMP)
  if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()

  set(MUMPS_LIBRARIES "")
  foreach(_lib dmumps mumps_common pord)
    find_library(MUMPS_${_lib}_LIBRARY NAMES ${_lib}
                 HINTS $ENV{MUMPS_DIR} PATH_SUFFIXES lib lib64)
    if(MUMPS_${_lib}_LIBRARY)
      list(APPEND MUMPS_LIBRARIES ${MUMPS_${_lib}_LIBRARY})
    elseif(NOT _lib STREQUAL "pord")   # pord optional
      message(FATAL_ERROR "MFEM built with MUMPS but lib${_lib} not found (set MUMPS_DIR)")
    endif()
  endforeach()
  message(STATUS "MUMPS libraries: ${MUMPS_LIBRARIES}")
endif(MFEM_USE_MUMPS)

# Petsc
if(MFEM_USE_PETSC)
  message(STATUS "MFEM: using PETSc, so SLOTH will use PETSc")
  find_package(PETSc REQUIRED)
endif(MFEM_USE_PETSC)

if(MFEM_USE_SUITESPARSE)
  find_package(SuiteSparse REQUIRED)
endif(MFEM_USE_SUITESPARSE)

find_package(METIS REQUIRED)

# HYPRE version
include(hypre)
set(MFEM_HYPRE_VERSION ${HYPRE_VERSION}) 

add_compile_definitions(MFEM_HYPRE_VERSION=${MFEM_HYPRE_VERSION})