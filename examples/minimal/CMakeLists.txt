# Minimal program using an installed SLOTH, organised as for a practical session:
# one sub-directory per case.
#   allen_cahn/  : 1D Allen-Cahn interface at equilibrium (tanh profile), whose
#                  double-well potential is generated from coefficients.json,
#                  checked against the analytical solution
# To add a case: copy allen_cahn/, rename it (and its targets/tests) and add an
# add_subdirectory() line at the end of this file.
#
# SLOTH must be visible to CMake:
#   - SLOTH installed with envSloth.sh --install=DIR + make install:
#         export CMAKE_PREFIX_PATH=DIR        (or cmake -DCMAKE_PREFIX_PATH=DIR)
#   - SLOTH installed with spack:
#         spack load sloth
#
# Build, run and check:
#   cmake -B build && cmake --build build
#   cd build && ctest --output-on-failure
# (or run the case directly: cd build/allen_cahn && mpirun -np 1 ./allen_cahn)
cmake_minimum_required(VERSION 3.13)
project(sloth_minimal LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

find_package(Sloth REQUIRED)

enable_testing()
find_program(MPIEXEC_EXECUTABLE NAMES mpiexec mpirun)

# Python of the active environment first (conda/mamba env, venv), then the
# first python3 found in PATH (as SLOTH itself does, see SYMPY4SLOTH.cmake).
# Needed with SymPy to generate the coefficients from the JSON files; the
# result checks only use the standard library.
set(Python3_FIND_VIRTUALENV FIRST)
set(Python3_FIND_STRATEGY LOCATION)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
message(STATUS "Python used for the coefficients and the checks: ${Python3_EXECUTABLE}")

execute_process(COMMAND ${Python3_EXECUTABLE} -c "import sympy"
                RESULT_VARIABLE _no_sympy OUTPUT_QUIET ERROR_QUIET)
if(NOT _no_sympy EQUAL 0)
  message(FATAL_ERROR "SymPy not found for ${Python3_EXECUTABLE}: install it (pip install sympy) "
                      "or activate an environment that provides it")
endif()
if(NOT EXISTS "${Sloth_COEFFICIENT_SCRIPT}")
  message(FATAL_ERROR "SLOTH coefficient generator not found ('${Sloth_COEFFICIENT_SCRIPT}')")
endif()

add_subdirectory(allen_cahn)
