# ##############################################################################
# SlothMFEM: INTERFACE target carrying everything needed to compile and link
# against MFEM, built from MFEM's own configuration file (config.mk):
#   - include directories : MFEM_INC_DIR + the -I paths of MFEM_TPLFLAGS
#   - link                : MFEM_LIB_FILE + MFEM_EXT_LIBS (HYPRE, METIS,
#                           SuiteSparse, MUMPS, PETSc, BLAS/LAPACK, ...,
#                           exactly as used to build MFEM) + MPI if needed
# The third-party libraries of MFEM no longer need to be searched one by one.
#
# Input : MFEM_CONFIG_FILE (path of config.mk, located by FindMFEM.cmake)
# Output: target SlothMFEM, variables SLOTHMFEM_INCLUDE_DIRS, SLOTHMFEM_LINK
# ##############################################################################

if(NOT EXISTS "${MFEM_CONFIG_FILE}")
  message(FATAL_ERROR "SlothMFEM: MFEM configuration file not found ('${MFEM_CONFIG_FILE}')")
endif()

# Re-run CMake automatically (at the next "make") when config.mk changes,
# e.g. after "brew upgrade mfem" or a re-installation of MFEM with spack
set_property(DIRECTORY "${PROJECT_SOURCE_DIR}" APPEND PROPERTY
             CMAKE_CONFIGURE_DEPENDS "${MFEM_CONFIG_FILE}")

# --- Read config.mk: "NAME = value" lines ------------------------------------
file(STRINGS "${MFEM_CONFIG_FILE}" _mk_lines REGEX "^[A-Za-z_][A-Za-z0-9_]*[ \t]*=")
set(_mk_names "")
foreach(_l IN LISTS _mk_lines)
  if(_l MATCHES "^([A-Za-z_][A-Za-z0-9_]*)[ \t]*=[ \t]*(.*)$")
    set(_mk_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}")
    list(APPEND _mk_names ${CMAKE_MATCH_1})
  endif()
endforeach()

# Expand make-style references $(NAME) using the values read above
function(_slothmfem_expand out value)
  foreach(_i RANGE 5)
    if(NOT value MATCHES "\\$\\(([A-Za-z_][A-Za-z0-9_]*)\\)")
      break()
    endif()
    string(REPLACE "$(${CMAKE_MATCH_1})" "${_mk_${CMAKE_MATCH_1}}" value "${value}")
  endforeach()
  string(STRIP "${value}" value)
  set(${out} "${value}" PARENT_SCOPE)
endfunction()

# Homebrew: replace version-specific paths ".../Cellar/<pkg>/<version>/..." by
# the stable ".../opt/<pkg>/..." when the latter exists and points to the same
# files, so that a "brew upgrade" does not break existing build directories.
# Paths outside a Homebrew Cellar are returned unchanged.
function(_slothmfem_stable_path out path)
  set(_res "${path}")
  if(path MATCHES "^(.*)/Cellar/([^/]+)/[^/]+(/.*)?$")
    set(_cand "${CMAKE_MATCH_1}/opt/${CMAKE_MATCH_2}${CMAKE_MATCH_3}")
    if(EXISTS "${_cand}")
      get_filename_component(_real_cand "${_cand}" REALPATH)
      get_filename_component(_real_path "${path}" REALPATH)
      if(_real_cand STREQUAL _real_path)
        set(_res "${_cand}")
      endif()
    endif()
  endif()
  set(${out} "${_res}" PARENT_SCOPE)
endfunction()

# Versioned library file name (e.g. libmfem.4.10.0.dylib, libmfem.so.4.10.0):
# use the unversioned name (libmfem.dylib / libmfem.so) of the same directory
# when it exists and points to the same file.
function(_slothmfem_unversioned_lib out file)
  set(_res "${file}")
  get_filename_component(_dir "${file}" DIRECTORY)
  get_filename_component(_name "${file}" NAME)
  if(_name MATCHES "^(lib[^.]+)\\.([0-9.]+\\.)?(dylib|so)(\\.[0-9.]+)?$")
    set(_cand "${_dir}/${CMAKE_MATCH_1}.${CMAKE_MATCH_3}")
    if(EXISTS "${_cand}")
      get_filename_component(_real_cand "${_cand}" REALPATH)
      get_filename_component(_real_file "${file}" REALPATH)
      if(_real_cand STREQUAL _real_file)
        set(_res "${_cand}")
      endif()
    endif()
  endif()
  set(${out} "${_res}" PARENT_SCOPE)
endfunction()

_slothmfem_expand(_mfem_inc_dir  "${_mk_MFEM_INC_DIR}")
_slothmfem_expand(_mfem_lib_file "${_mk_MFEM_LIB_FILE}")
_slothmfem_expand(_mfem_tplflags "${_mk_MFEM_TPLFLAGS}")
_slothmfem_expand(_mfem_ext_libs "${_mk_MFEM_EXT_LIBS}")

if(NOT EXISTS "${_mfem_lib_file}")
  message(FATAL_ERROR "SlothMFEM: MFEM library '${_mfem_lib_file}' (MFEM_LIB_FILE in ${MFEM_CONFIG_FILE}) does not exist")
endif()
_slothmfem_stable_path(_mfem_inc_dir "${_mfem_inc_dir}")
_slothmfem_stable_path(_mfem_lib_file "${_mfem_lib_file}")
_slothmfem_unversioned_lib(_mfem_lib_file "${_mfem_lib_file}")

# --- Include directories -----------------------------------------------------
set(SLOTHMFEM_INCLUDE_DIRS "${_mfem_inc_dir}")
separate_arguments(_tpl UNIX_COMMAND "${_mfem_tplflags}")
set(_next_is_dir FALSE)
foreach(_f IN LISTS _tpl)
  if(_next_is_dir)
    list(APPEND SLOTHMFEM_INCLUDE_DIRS "${_f}")
    set(_next_is_dir FALSE)
  elseif(_f STREQUAL "-I" OR _f STREQUAL "-isystem")
    set(_next_is_dir TRUE)
  elseif(_f MATCHES "^-I(.+)$")
    list(APPEND SLOTHMFEM_INCLUDE_DIRS "${CMAKE_MATCH_1}")
  elseif(_f MATCHES "^-isystem(.+)$")
    # (separate branch: with "A MATCHES ... OR B MATCHES ...", the second test
    #  resets CMAKE_MATCH_1)
    list(APPEND SLOTHMFEM_INCLUDE_DIRS "${CMAKE_MATCH_1}")
  endif()
endforeach()
set(_stable_incs "")
foreach(_d IN LISTS SLOTHMFEM_INCLUDE_DIRS)
  _slothmfem_stable_path(_d "${_d}")
  list(APPEND _stable_incs "${_d}")
endforeach()
set(SLOTHMFEM_INCLUDE_DIRS ${_stable_incs})
list(REMOVE_DUPLICATES SLOTHMFEM_INCLUDE_DIRS)

# --- Link items --------------------------------------------------------------
# MFEM_EXT_LIBS is kept in order (order matters for static libraries); only
# repeated -L<dir> and -Wl,-rpath,<dir> entries are dropped.
separate_arguments(_ext UNIX_COMMAND "${_mfem_ext_libs}")
set(SLOTHMFEM_LINK "${_mfem_lib_file}")
set(_seen_dirs "")
foreach(_f IN LISTS _ext)
  # Stable Homebrew paths for -L<dir>, -Wl,-rpath,<dir> and full library paths
  if(_f MATCHES "^-L(.+)$")
    _slothmfem_stable_path(_d "${CMAKE_MATCH_1}")
    set(_f "-L${_d}")
  elseif(_f MATCHES "^-Wl,-rpath,(.+)$")
    _slothmfem_stable_path(_d "${CMAKE_MATCH_1}")
    set(_f "-Wl,-rpath,${_d}")
  elseif(_f MATCHES "^/")
    _slothmfem_stable_path(_f "${_f}")
  endif()
  if(_f MATCHES "^-L" OR _f MATCHES "^-Wl,-rpath,")
    if(_f IN_LIST _seen_dirs)
      continue()
    endif()
    list(APPEND _seen_dirs "${_f}")
  endif()
  list(APPEND SLOTHMFEM_LINK "${_f}")
endforeach()

# --- Target ------------------------------------------------------------------
if(NOT TARGET SlothMFEM)
  add_library(SlothMFEM INTERFACE)
endif()
target_include_directories(SlothMFEM SYSTEM INTERFACE ${SLOTHMFEM_INCLUDE_DIRS})
target_link_libraries(SlothMFEM INTERFACE ${SLOTHMFEM_LINK})
if(MFEM_USE_MPI)
  # MPI may not appear in MFEM_EXT_LIBS when MFEM was built with mpicxx.
  # Plain paths (not the MPI::MPI_CXX target), so that the exported
  # SlothTargets.cmake stays self-contained (no find_package(MPI) needed by
  # projects using an installed Sloth), as with the previous link line.
  set(SLOTHMFEM_MPI_INCLUDE_DIRS "")
  foreach(_d IN LISTS MPI_CXX_INCLUDE_DIRS)
    _slothmfem_stable_path(_d "${_d}")
    list(APPEND SLOTHMFEM_MPI_INCLUDE_DIRS "${_d}")
  endforeach()
  set(SLOTHMFEM_MPI_LIBRARIES "")
  foreach(_l IN LISTS MPI_CXX_LIBRARIES)
    _slothmfem_stable_path(_l "${_l}")
    list(APPEND SLOTHMFEM_MPI_LIBRARIES "${_l}")
  endforeach()
  target_include_directories(SlothMFEM SYSTEM INTERFACE ${SLOTHMFEM_MPI_INCLUDE_DIRS})
  target_link_libraries(SlothMFEM INTERFACE ${SLOTHMFEM_MPI_LIBRARIES})
endif()

message(STATUS "SlothMFEM: MFEM library ${_mfem_lib_file}")
