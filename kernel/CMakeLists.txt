# Create the library (SLOTH_* options are declared in the top-level CMakeLists.txt)
if(SLOTH_USE_SHARED)
    add_library(Sloth SHARED)
else()
    add_library(Sloth STATIC)
endif(SLOTH_USE_SHARED)

# Get all cpp files
file(GLOB_RECURSE SLOTH_SRC CONFIGURE_DEPENDS
    ${CMAKE_SOURCE_DIR}/kernel/src/*.cpp
)

# Exclude specific cpp files if not required 
if(NOT SLOTH_USE_EXPRTK)
    list(REMOVE_ITEM SLOTH_SRC ${CMAKE_SOURCE_DIR}/kernel/src/Coefficients/ExprTkSloth.cpp)
endif()

target_sources(Sloth PRIVATE ${SLOTH_SRC})
set_compile_options(Sloth)

set_target_properties(Sloth PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CMAKE_CXX_EXTENSIONS OFF
)

# And the most important part set the include directories
target_include_directories(Sloth
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/kernel/include>
        $<INSTALL_INTERFACE:include>
)

# ExprTk
# ExprTk (header-only): directory containing exprtk.hpp, stored in the cache
# variable ExprTk_INCLUDE_DIR.
#   - envSloth.sh --exprtk=DIR exports the environment variable
#     ExprTk_INCLUDE_DIR, which takes priority (a new --exprtk path in an
#     existing build directory is taken into account);
#   - otherwise the cached value is kept (a plain "cmake ." keeps ExprTk);
#   - -DExprTk_INCLUDE_DIR=DIR can also be given directly.
if(SLOTH_USE_EXPRTK)
    if(NOT "$ENV{ExprTk_INCLUDE_DIR}" STREQUAL ""
       AND NOT "$ENV{ExprTk_INCLUDE_DIR}" STREQUAL "${ExprTk_INCLUDE_DIR}")
        unset(ExprTk_INCLUDE_DIR CACHE)
    endif()
    find_path(ExprTk_INCLUDE_DIR
        NAMES exprtk.hpp
        HINTS ENV ExprTk_INCLUDE_DIR
        PATH_SUFFIXES exprtk include/exprtk
        DOC "Directory containing exprtk.hpp"
    )
    if(NOT ExprTk_INCLUDE_DIR)
        message(FATAL_ERROR
            "SLOTH_USE_EXPRTK is ON but exprtk.hpp was not found.\n"
            "Use envSloth.sh --exprtk=<dir containing exprtk.hpp> "
            "or -DExprTk_INCLUDE_DIR=<dir>.")
    endif()
    target_include_directories(Sloth SYSTEM PUBLIC
        $<BUILD_INTERFACE:${ExprTk_INCLUDE_DIR}>
    )
endif()


# MFEM and all its third-party libraries (from MFEM's config.mk, see
# cmake/modules/SLOTHMFEM4SLOTH.cmake)
target_link_libraries(Sloth PUBLIC SlothMFEM)


# ########################
# MUMPS ??
# ########################
# (MUMPS libraries are linked through SlothMFEM)
if(MFEM_USE_MUMPS)
  target_compile_definitions(Sloth PUBLIC SLOTH_USE_MUMPS=1)
endif()

# EXTERNAL??
# External library (envSloth.sh --external=EXT_LIBDIR,EXT_LIBNAME)
#   SLOTH_EXTERNAL_LIBDIR  : directory containing the library
#   SLOTH_EXTERNAL_LIBNAME : library name, as given to the linker
# Linking is unchanged: link directory + library name.
# The values are stored in the cache: the environment variables EXT_LIBDIR /
# EXT_LIBNAME (exported by envSloth.sh) take priority; otherwise the cached
# values are kept (a plain "cmake ." keeps the external library). They can also
# be given directly with -DSLOTH_EXTERNAL_LIBDIR=... -DSLOTH_EXTERNAL_LIBNAME=...
if(SLOTH_USE_EXTERNAL)
    set(SLOTH_EXTERNAL_LIBDIR "" CACHE PATH "External library: directory (EXT_LIBDIR)")
    set(SLOTH_EXTERNAL_LIBNAME "" CACHE STRING "External library: name given to the linker (EXT_LIBNAME)")
    if(NOT "$ENV{EXT_LIBDIR}" STREQUAL "")
        set(SLOTH_EXTERNAL_LIBDIR "$ENV{EXT_LIBDIR}" CACHE PATH "External library: directory (EXT_LIBDIR)" FORCE)
    endif()
    if(NOT "$ENV{EXT_LIBNAME}" STREQUAL "")
        set(SLOTH_EXTERNAL_LIBNAME "$ENV{EXT_LIBNAME}" CACHE STRING "External library: name given to the linker (EXT_LIBNAME)" FORCE)
    endif()

    if(SLOTH_EXTERNAL_LIBDIR STREQUAL "" OR NOT IS_DIRECTORY "${SLOTH_EXTERNAL_LIBDIR}")
        message(FATAL_ERROR
            "SLOTH_USE_EXTERNAL is ON but the external library directory "
            "'${SLOTH_EXTERNAL_LIBDIR}' does not exist.\n"
            "Use envSloth.sh --external=EXT_LIBDIR,EXT_LIBNAME "
            "or -DSLOTH_EXTERNAL_LIBDIR=<dir>.")
    endif()
    if(SLOTH_EXTERNAL_LIBNAME STREQUAL "")
        message(FATAL_ERROR
            "SLOTH_USE_EXTERNAL is ON but no external library name was given.\n"
            "Use envSloth.sh --external=EXT_LIBDIR,EXT_LIBNAME "
            "or -DSLOTH_EXTERNAL_LIBNAME=<name>.")
    endif()

    # Early check (warning only, linking is left unchanged): does the library
    # exist in SLOTH_EXTERNAL_LIBDIR? Accepts "name", "libname.*" or a file name.
    file(GLOB _sloth_ext_candidates
         "${SLOTH_EXTERNAL_LIBDIR}/${SLOTH_EXTERNAL_LIBNAME}"
         "${SLOTH_EXTERNAL_LIBDIR}/lib${SLOTH_EXTERNAL_LIBNAME}.*")
    if(NOT _sloth_ext_candidates AND NOT EXISTS "${SLOTH_EXTERNAL_LIBNAME}")
        message(WARNING
            "External library '${SLOTH_EXTERNAL_LIBNAME}' not found in "
            "'${SLOTH_EXTERNAL_LIBDIR}' (looked for '${SLOTH_EXTERNAL_LIBNAME}' "
            "and 'lib${SLOTH_EXTERNAL_LIBNAME}.*'). Linking may fail.")
    endif()

    target_link_directories(Sloth PUBLIC ${SLOTH_EXTERNAL_LIBDIR})
    target_link_libraries(Sloth PUBLIC ${SLOTH_EXTERNAL_LIBNAME})
endif(SLOTH_USE_EXTERNAL)


# LibTorch API ??
if(SLOTH_USE_LIBTORCH)
    find_package(Torch REQUIRED)
    MESSAGE(STATUS "Sloth use libTorch API ")
    target_include_directories(Sloth PUBLIC ${TORCH_INCLUDE_DIRS})

    target_link_libraries(Sloth PUBLIC ${TORCH_LIBRARIES})
    
    target_compile_definitions(Sloth PUBLIC SLOTH_USE_LIBTORCH=1)
endif(SLOTH_USE_LIBTORCH)

# ########################
# VORO++ API ??
# ########################
if(SLOTH_USE_VORO)
    # Search order:
    #   1. $ENV{VORO} (set by envSloth.sh --voro=DIR): DIR/include, DIR/lib
    #   2. Homebrew / /usr/local
    #   3. CMake default locations (CMAKE_PREFIX_PATH, system)
    # Results are cached (VORO_INCLUDE_DIR, VORO_LIBRARY): a plain "cmake ."
    # keeps them. When --voro points to another directory than at the
    # previous configure, the search is done again.
    # (Previously $ENV{VORO} came after the Homebrew paths, so --voro was
    # ignored whenever voro++ was also installed with Homebrew.)
    if(NOT "$ENV{VORO}" STREQUAL "" AND NOT "$ENV{VORO}" STREQUAL "${SLOTH_VORO_ROOT_USED}")
        unset(VORO_INCLUDE_DIR CACHE)
        unset(VORO_LIBRARY CACHE)
        set(SLOTH_VORO_ROOT_USED "$ENV{VORO}" CACHE INTERNAL "Voro++ root (--voro) used for the cached search")
    endif()

    set(_voro_inc_hints /opt/homebrew/include /usr/local/include)
    set(_voro_lib_hints /opt/homebrew/lib /usr/local/lib)
    if(NOT "$ENV{VORO}" STREQUAL "")
        list(INSERT _voro_inc_hints 0 "$ENV{VORO}/include")
        list(INSERT _voro_lib_hints 0 "$ENV{VORO}/lib")
    endif()

    find_path(
        VORO_INCLUDE_DIR
        NAMES
        voro++.hh
        voro++/voro++.hh
        HINTS ${_voro_inc_hints}
        DOC "Voro++ include directory"
    )

    find_library(
        VORO_LIBRARY
        NAMES voro++ libvoro++
        HINTS ${_voro_lib_hints}
        DOC "Voro++ library"
    )

    message(STATUS "VORO_INCLUDE_DIR=${VORO_INCLUDE_DIR}")
    message(STATUS "VORO_LIBRARY=${VORO_LIBRARY}")
    if(NOT VORO_INCLUDE_DIR)
        message(FATAL_ERROR "Voro++ include not found")
    endif()

    if(EXISTS "${VORO_INCLUDE_DIR}/voro++.hh")
        target_include_directories(Sloth PUBLIC ${VORO_INCLUDE_DIR})
    elseif(EXISTS "${VORO_INCLUDE_DIR}/voro++/voro++.hh")
        target_include_directories(Sloth PUBLIC ${VORO_INCLUDE_DIR}/voro++)
    endif()

    if(NOT VORO_LIBRARY)
        message(FATAL_ERROR "Voro++ library not found")
    endif()

    target_link_libraries(Sloth PUBLIC ${VORO_LIBRARY})
    
    message(STATUS "Sloth use Voro++ library")
endif(SLOTH_USE_VORO)



# ExprTk ??
if(SLOTH_USE_EXPRTK)
    target_compile_definitions(Sloth PUBLIC SLOTH_USE_EXPRTK=1)
endif()

## Debug option
if(MATOOLS_VERBOSITY_LEVEL_1)
  target_compile_definitions(Sloth PUBLIC MATIMERS_VEROBSITY_LEVEL_1)
endif()


# ########################
# # Installation 
# ########################
install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION include
)

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/scripts/
    DESTINATION scripts
)
# SlothTargets.cmake creation
# SlothMFEM: exported with Sloth, which depends on it
install(TARGETS SlothMFEM EXPORT SlothTargets)

install(
    TARGETS Sloth
    EXPORT SlothTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include
)



install(
    EXPORT SlothTargets
    FILE SlothTargets.cmake
    NAMESPACE Sloth::
    DESTINATION cmake/Sloth
)

