cmake_minimum_required(VERSION 3.12)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON )

project(sloth VERSION 2.1 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ########################
# # SLOTH options
# ########################
# All options are declared here, before being used, so that they are listed
# (with their description) by "cmake -LH" or "ccmake".
# envSloth.sh sets them through -D... (see its command-line options).
option(SLOTH_USE_SHARED   "Build Sloth as a shared library (default: static)" OFF)
option(SLOTH_USE_EXPRTK   "Enable ExprTk (analytical expressions parsed at runtime)" OFF)
option(SLOTH_USE_VORO     "Enable Voro++ (Voronoi tessellations)" OFF)
option(SLOTH_USE_LIBTORCH "Enable the LibTorch API" OFF)
option(SLOTH_USE_EXTERNAL "Link Sloth with an external library (EXT_LIBDIR, EXT_LIBNAME)" OFF)
option(MATOOLS_VERBOSITY_LEVEL_1 "MATools profiling: verbosity level 1" OFF)
# Optimisation for the build machine (Release / RelWithDebInfo / MinSizeRel).
# Set them to OFF when SLOTH is built for other machines or used with another
# compiler (spack package, shared installation): -march=native code may not run
# on another CPU, and LTO objects can only be linked with the same compiler.
option(SLOTH_NATIVE_ARCH "Optimise for the CPU of the build machine (-march=native)" ON)
option(SLOTH_ENABLE_LTO  "Link-time optimisation (Release, MinSizeRel)" ${SLOTH_NATIVE_ARCH})

## Force Release Mode
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
             Release PortableRelease RelWithDebInfo MinSizeRel Optim Debug Coverage)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()



# ########################
# # Global Configuration
# ########################
# SlothInstall : default name of the repository for installation
# Priority:
#   1. USER_INSTALL_PATH environment variable (set by envSloth.sh --install=...)
#   2. -DCMAKE_INSTALL_PREFIX=... given explicitly by the user
#   3. default: ../SlothInstall
# The value is stored in the cache, so that the cache, the configuration
# summary and "make install" all use the same directory.
# The default (3) is applied:
#   - at the first configure (or "/usr/local", CMake's own default, found in
#     the cache of build directories created before this change);
#   - when envSloth.sh is run again without --install (it always exports
#     USER_INSTALL_PATH, empty in that case) and the cached prefix was set
#     automatically (by a previous --install or by the default), as before.
# A plain "cmake ." (USER_INSTALL_PATH not in the environment at all) keeps
# the cached prefix, whatever it is.
if(NOT "$ENV{USER_INSTALL_PATH}" STREQUAL "")
    set(_sloth_auto_prefix "$ENV{USER_INSTALL_PATH}")
elseif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR CMAKE_INSTALL_PREFIX STREQUAL "/usr/local")
    set(_sloth_auto_prefix "${CMAKE_SOURCE_DIR}/../SlothInstall")
elseif(DEFINED ENV{USER_INSTALL_PATH}
       AND CMAKE_INSTALL_PREFIX STREQUAL "${SLOTH_AUTO_INSTALL_PREFIX}")
    set(_sloth_auto_prefix "${CMAKE_SOURCE_DIR}/../SlothInstall")
else()
    set(_sloth_auto_prefix "")
endif()

if(NOT _sloth_auto_prefix STREQUAL "")
    set(CMAKE_INSTALL_PREFIX "${_sloth_auto_prefix}" CACHE PATH "Installation directory" FORCE)
    set(SLOTH_AUTO_INSTALL_PREFIX "${_sloth_auto_prefix}" CACHE INTERNAL "Install prefix set automatically by SLOTH")
endif()

# ########################
# # SymPy 
# ########################
include(SYMPY4SLOTH)

check_python_module(sympy SYMPY_FOUND)
if(NOT SYMPY_FOUND)
    message(FATAL_ERROR "SymPy is required to build Sloth. Please install SymPy before building Sloth.")
endif()


# ########################
# # MFEM Configuration
# ########################
include(MFEM4SLOTH)

# ########################
# # Compiler options
# ########################
include(COMPILE4SLOTH)

# ########################
# # Sloth kernel
# ########################
add_subdirectory(kernel)

# ########################
# # Doc Configuration
# ########################
if(NOT TARGET doc)
  add_subdirectory(doc)
endif(NOT TARGET doc)


# ########################
# LibTorch API ??
# ########################
if(SLOTH_USE_LIBTORCH)
    find_package(Torch REQUIRED)
endif(SLOTH_USE_LIBTORCH)



# ########################
# Script
# ########################
set(SCRIPT_FILE ${CMAKE_SOURCE_DIR}/scripts/GenerateCoefficient.py)
set(BUILD_SCRIPT_DIR ${CMAKE_BINARY_DIR}/scripts)
file(COPY ${SCRIPT_FILE} DESTINATION ${BUILD_SCRIPT_DIR})

# ########################
# Tests
# ########################
# By default, make command without target only builds Sloth
add_custom_target(sloth COMMENT "Build Sloth " DEPENDS Sloth)

include(TEST4SLOTH)

include(CTest)
enable_testing()
add_subdirectory(tests EXCLUDE_FROM_ALL)

# By default, build all tests only after Sloth has been compiled.
# Otherwise, Sloth is compiled and then, all tests are built.
add_custom_target(sloth_tests COMMENT "Build Sloth + all tests " DEPENDS Sloth tests)



# ########################
# # Code Coverage
# ########################
IF(CMAKE_BUILD_TYPE MATCHES Coverage)
  # Collect CodeCoverage analysis with gcovr
  add_custom_target(cc
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/Coverage
    COMMAND gcovr -r ${CMAKE_BINARY_DIR}/kernel --filter ${CMAKE_SOURCE_DIR}/kernel/src --filter ${CMAKE_SOURCE_DIR}/kernel/include  --exclude '${CMAKE_SOURCE_DIR}/kernel/include/MAToolsProfiling/.*'  --exclude-unreachable-branches --xml=${CMAKE_BINARY_DIR}/Coverage/coverage.xml --html --html-details -o ${CMAKE_BINARY_DIR}/Coverage/coverage.html
    COMMENT "Generate Html report for code coverage analysis (see Coverage/coverage.html)"
  )
ENDIF()

# ########################
# # Static analysis
# ########################
# Check if the lint target is already defined
if(NOT TARGET lint)
  find_program(CPP_LINT cpplint)

  # Add a custom target for linting
  add_custom_target(lint
    COMMAND ${CPP_LINT}
    --filter=-runtime/references,-build/header_guard,-runtime/string,-build/include
    --counting=total
    --extensions=cpp,cxx,hpp,h,hxx,tpp
    --linelength=100
    --exclude ${CMAKE_SOURCE_DIR}/kernel/MAToolsProfiling
    --recursive ${CMAKE_SOURCE_DIR}/kernel ${CMAKE_SOURCE_DIR}/tests
    COMMENT "Linting source files"
  )
endif()

# ########################
# # Configuration summary (read-only)
# ########################
include(SUMMARY4SLOTH)
sloth_print_summary()
