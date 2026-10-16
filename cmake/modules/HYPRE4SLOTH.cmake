# ##############################################################################
# HYPRE version
#
# Sets the cache variable HYPRE_VERSION (integer, e.g. 30200 for HYPRE 3.2.0),
# passed to the compiler as MFEM_HYPRE_VERSION (see MFEM4SLOTH.cmake).
#
# 1. Read from MFEM's own configuration header (MFEM_HYPRE_VERSION in
#    mfem/config/_config.hpp): MFEM was built against this HYPRE version.
# 2. Fallback: read HYPRE_RELEASE_VERSION (or HYPRE_PACKAGE_VERSION) from
#    HYPRE_config.h, looked for in the include directories of MFEM's
#    config.mk (SLOTHMFEM_INCLUDE_DIRS), then in $ENV{HYPRE_DIR}/include.
#    Conversion: "X.Y.Z" -> X*10000 + Y*100 + Z (e.g. "2.10.0b" -> 21000).
#
# No program is compiled: HYPRE_DIR is not required when MFEM provides the
# version, and a plain "cmake ." in an existing build directory works.
# ##############################################################################

if(MFEM_USE_MPI)
  unset(_sloth_hypre_version)

  # 1. MFEM configuration header
  foreach(_mfem_inc IN LISTS MFEM_INCLUDE_DIRS)
    foreach(_cfg mfem/config/_config.hpp config/_config.hpp
                 mfem/config/config.hpp  config/config.hpp)
      if(NOT DEFINED _sloth_hypre_version AND EXISTS "${_mfem_inc}/${_cfg}")
        file(STRINGS "${_mfem_inc}/${_cfg}" _line
             REGEX "^[ \t]*#[ \t]*define[ \t]+MFEM_HYPRE_VERSION[ \t]+[0-9]+")
        if(_line)
          string(REGEX REPLACE ".*MFEM_HYPRE_VERSION[ \t]+([0-9]+).*" "\\1"
                 _sloth_hypre_version "${_line}")
          message(STATUS "HYPRE version read from ${_mfem_inc}/${_cfg}")
        endif()
      endif()
    endforeach()
  endforeach()

  # 2. Fallback: HYPRE_config.h
  if(NOT DEFINED _sloth_hypre_version)
    set(_hypre_incs ${SLOTHMFEM_INCLUDE_DIRS})
    if(NOT "$ENV{HYPRE_DIR}" STREQUAL "")
      list(APPEND _hypre_incs "$ENV{HYPRE_DIR}/include")
    endif()
    foreach(_hypre_inc IN LISTS _hypre_incs)
      if(NOT DEFINED _sloth_hypre_version AND EXISTS "${_hypre_inc}/HYPRE_config.h")
        file(STRINGS "${_hypre_inc}/HYPRE_config.h" _line
             REGEX "^[ \t]*#[ \t]*define[ \t]+HYPRE_(RELEASE|PACKAGE)_VERSION[ \t]")
        # First 2 or 3 numeric components: "3.2.0", "2.10.0b", "2.31"
        if(_line MATCHES "\"([0-9]+)\\.([0-9]+)(\\.([0-9]+))?")
          set(_maj ${CMAKE_MATCH_1})
          set(_min ${CMAKE_MATCH_2})
          set(_pat ${CMAKE_MATCH_4})
          if("${_pat}" STREQUAL "")
            set(_pat 0)
          endif()
          math(EXPR _sloth_hypre_version "${_maj}*10000 + ${_min}*100 + ${_pat}")
          message(STATUS "HYPRE version read from ${_hypre_inc}/HYPRE_config.h")
        endif()
      endif()
    endforeach()
  endif()

  if(NOT DEFINED _sloth_hypre_version)
    message(FATAL_ERROR
      "Unable to determine the HYPRE version: MFEM_HYPRE_VERSION is not defined "
      "in the MFEM configuration header and HYPRE_config.h was not found "
      "(searched in: ${_hypre_incs}). Set the HYPRE_DIR environment variable.")
  endif()

  set(HYPRE_VERSION ${_sloth_hypre_version} CACHE STRING "HYPRE version." FORCE)
  message(STATUS "Found HYPRE version ${HYPRE_VERSION}")
endif(MFEM_USE_MPI)
