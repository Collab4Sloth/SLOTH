cmake_minimum_required(VERSION 3.12)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON )

project(sloth VERSION 1.0 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

## Force Release Mode
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Coverage RelWithDebInfo Debug)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# Check GCC version for including filesystem
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    find_package(Filesystem REQUIRED)
    set(FILESYSTEM_LIB stdc++fs)
  endif()
endif()


# ########################
# # Global Configuration
# ########################
# SlothInstall : default name of the repository for installation
if("$ENV{USER_INSTALL_PATH}" STREQUAL "")
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/../SlothInstall" )
else()
    set(CMAKE_INSTALL_PREFIX "$ENV{USER_INSTALL_PATH}")
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
# # Test Configuration
# ########################
include(COMPILE4SLOTH)

# ########################
# # Sloth header's only
# # library (interface)
# ########################
# include_directories(${CMAKE_SOURCE_DIR})
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
option(SLOTH_USE_LIBTORCH "Enable the use of libtorch API" OFF)
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
