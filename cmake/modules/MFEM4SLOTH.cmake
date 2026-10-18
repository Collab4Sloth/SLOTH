# ##############################################################################
# MFEM and its third-party libraries: everything comes from MFEM's own
# configuration (config.mk), see SLOTHMFEM4SLOTH.cmake (target SlothMFEM).
# ##############################################################################
include(SLOTHMFEM4SLOTH)

if(MFEM_USE_MPI)
  message(STATUS "MFEM: using MPI, so SLOTH will use MPI (with HYPRE, METIS)")
endif(MFEM_USE_MPI)

# Support for OpenMP
if(MFEM_USE_OPENMP OR MFEM_USE_LEGACY_OPENMP)
  message(STATUS "MFEM: using OPENMP, so SLOTH will use OPENMP")
  find_package(OpenMP REQUIRED)

  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    # Exported with SlothMFEM so that applications using find_package(Sloth) link OpenMP too
    separate_arguments(_sloth_omp_flags UNIX_COMMAND "${OpenMP_CXX_FLAGS}")
    target_compile_options(SlothMFEM INTERFACE ${_sloth_omp_flags})
    target_link_libraries(SlothMFEM INTERFACE ${OpenMP_CXX_LIBRARIES})
  endif()
endif(MFEM_USE_OPENMP OR MFEM_USE_LEGACY_OPENMP)

# MUMPS: libraries come from MFEM_EXT_LIBS (SlothMFEM)
if(MFEM_USE_MUMPS)
  message(STATUS "MFEM: using MUMPS, so SLOTH will use MUMPS")

  # MUMPS (spack) may be built with OpenMP
  find_package(OpenMP)
  if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    # Exported with SlothMFEM so that applications using find_package(Sloth) link the OpenMP runtime
    target_link_libraries(SlothMFEM INTERFACE ${OpenMP_CXX_LIBRARIES})
  endif()
endif(MFEM_USE_MUMPS)

# HYPRE version
include(HYPRE4SLOTH)
set(MFEM_HYPRE_VERSION ${HYPRE_VERSION})

add_compile_definitions(MFEM_HYPRE_VERSION=${MFEM_HYPRE_VERSION})
