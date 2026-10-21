# ##############################################################################
# MFEM for SLOTH, from MFEM's own configuration file (config.mk)
#
# 1. Locate config.mk (same search as the former FindMFEM.cmake):
#    CMAKE_PREFIX_PATH (environment: spack load / spack env / module load),
#    then MFEM_DIR (CMake or environment variable, set by envSloth.sh:
#    spack location -i mfem, brew --prefix mfem, or --mfem=DIR),
#    looking for <dir>/share/mfem/config.mk or <dir>/config.mk.
# 2. Read it: MFEM_USE_* (YES/NO), include directories, library, and the
#    third-party libraries exactly as used to build MFEM (MFEM_EXT_LIBS:
#    HYPRE, METIS, MUMPS, PETSc, SUNDIALS, BLAS/LAPACK, ...).
# 3. Build the INTERFACE target SlothMFEM, linked by Sloth.
#
# Outputs: target SlothMFEM, MFEM_CONFIG_FILE, MFEM_USE_*, MFEM_INCLUDE_DIRS,
#          MFEM_LIBRARY, SLOTHMFEM_INCLUDE_DIRS, SLOTHMFEM_LINK,
#          SLOTHMFEM_MPI_INCLUDE_DIRS, SLOTHMFEM_MPI_LIBRARIES
# ##############################################################################

# --- 1. Locate config.mk -----------------------------------------------------
if(MFEM_CONFIG_FILE AND NOT EXISTS "${MFEM_CONFIG_FILE}")
  # cached path no longer valid (MFEM removed or re-installed elsewhere)
  unset(MFEM_CONFIG_FILE CACHE)
endif()
string(REPLACE ":" ";" _sloth_env_prefix_path "$ENV{CMAKE_PREFIX_PATH}")
find_file(MFEM_CONFIG_FILE config.mk
  PATHS ${_sloth_env_prefix_path}
        ${MFEM_DIR}/build $ENV{MFEM_DIR}/build ${MFEM_DIR} $ENV{MFEM_DIR}
  PATH_SUFFIXES share/mfem
  NO_DEFAULT_PATH
  DOC "MFEM configuration file (config.mk)")
if(NOT MFEM_CONFIG_FILE)
  message(FATAL_ERROR
    "MFEM configuration file (share/mfem/config.mk) not found.\n"
    "Make MFEM visible with 'spack load mfem' / a spack environment / "
    "'module load', or use envSloth.sh --mfem=DIR, or -DMFEM_DIR=<MFEM prefix>.")
endif()
message(STATUS "Found MFEM: ${MFEM_CONFIG_FILE}")

# Re-run CMake automatically (at the next "make") when config.mk changes,
# e.g. after "brew upgrade mfem" or a re-installation of MFEM with spack
set_property(DIRECTORY "${PROJECT_SOURCE_DIR}" APPEND PROPERTY
             CMAKE_CONFIGURE_DEPENDS "${MFEM_CONFIG_FILE}")

# --- 2. Read config.mk: "NAME = value" lines ---------------------------------
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
# the stable ".../opt/<pkg>/..." when the latter exists and either points to the
# same files, or the versioned path no longer exists (MFEM bottle built against
# another version of the dependency, or dependency upgraded since), so that a
# "brew upgrade" does not break existing build directories or installed SLOTH.
# Paths outside a Homebrew Cellar are returned unchanged.
function(_slothmfem_stable_path out path)
  set(_res "${path}")
  if(path MATCHES "^(.*)/Cellar/([^/]+)/[^/]+(/.*)?$")
    set(_cand "${CMAKE_MATCH_1}/opt/${CMAKE_MATCH_2}${CMAKE_MATCH_3}")
    if(EXISTS "${_cand}")
      if(NOT EXISTS "${path}")
        set(_res "${_cand}")
      else()
        get_filename_component(_real_cand "${_cand}" REALPATH)
        get_filename_component(_real_path "${path}" REALPATH)
        if(_real_cand STREQUAL _real_path)
          set(_res "${_cand}")
        endif()
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

# MFEM_USE_XXX variables (YES / NO) for the rest of the configuration
foreach(_n IN LISTS _mk_names)
  if(_n MATCHES "^MFEM_USE_")
    _slothmfem_expand(${_n} "${_mk_${_n}}")
  endif()
endforeach()


# MPI (MFEM parallel build)
if(MFEM_USE_MPI)
  message(STATUS "MFEM uses MPI")
  find_package(MPI REQUIRED)
endif()

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
set(MFEM_INCLUDE_DIRS "${_mfem_inc_dir}")
set(MFEM_LIBRARY "${_mfem_lib_file}")

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
# Non-existent directories (stale paths recorded in config.mk) are dropped: they
# are useless for compiling and CMake rejects them in imported targets
# (find_package(Sloth) in applications)
set(_stable_incs "")
foreach(_d IN LISTS SLOTHMFEM_INCLUDE_DIRS)
  _slothmfem_stable_path(_d "${_d}")
  if(IS_DIRECTORY "${_d}")
    list(APPEND _stable_incs "${_d}")
  else()
    message(STATUS "SlothMFEM: include directory '${_d}' (from MFEM config.mk) does not exist, ignored")
  endif()
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

# --- 3. Target ---------------------------------------------------------------
if(NOT TARGET SlothMFEM)
  add_library(SlothMFEM INTERFACE)
endif()
target_include_directories(SlothMFEM SYSTEM INTERFACE ${SLOTHMFEM_INCLUDE_DIRS})
target_link_libraries(SlothMFEM INTERFACE ${SLOTHMFEM_LINK})
if(MFEM_USE_MPI)
  # MPI may not appear in MFEM_EXT_LIBS when MFEM was built with mpicxx.
  # Plain paths (not the MPI::MPI_CXX target), so that the exported
  # SlothTargets.cmake stays self-contained (no find_package(MPI) needed by
  # projects using an installed Sloth).
  set(SLOTHMFEM_MPI_INCLUDE_DIRS "")
  foreach(_d IN LISTS MPI_CXX_INCLUDE_DIRS)
    _slothmfem_stable_path(_d "${_d}")
    if(IS_DIRECTORY "${_d}")
      list(APPEND SLOTHMFEM_MPI_INCLUDE_DIRS "${_d}")
    endif()
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
