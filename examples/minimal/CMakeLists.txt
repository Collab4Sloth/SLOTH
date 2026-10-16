# Minimal program using an installed SLOTH.
#
# SLOTH must be visible to CMake:
#   - SLOTH installed with envSloth.sh --install=DIR + make install:
#         export CMAKE_PREFIX_PATH=DIR        (or cmake -DCMAKE_PREFIX_PATH=DIR)
#   - SLOTH installed with spack:
#         spack load sloth
#
# Build, run and check against the reference results:
#   cmake -B build && cmake --build build
#   cd build && ctest --output-on-failure
# (or run directly: mpirun -np 1 ./build/main)
cmake_minimum_required(VERSION 3.13)
project(sloth_minimal LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

find_package(Sloth REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE Sloth::Sloth)

# ------------------------------------------------------------------------------
# Tests (ctest): run the program, then compare its results with the reference
# results (ref/time_specialized.csv). The comparison needs Python 3 (standard
# library only); without it, only the run test is defined.
# ------------------------------------------------------------------------------
enable_testing()
find_program(MPIEXEC_EXECUTABLE NAMES mpiexec mpirun)
add_test(NAME minimal_run COMMAND ${MPIEXEC_EXECUTABLE} -n 1 $<TARGET_FILE:main>)
set_tests_properties(minimal_run PROPERTIES FIXTURES_SETUP minimal_results)

# Python of the active environment first (conda/mamba env, venv), then the
# first python3 found in PATH (as SLOTH itself does, see SYMPY4SLOTH.cmake)
set(Python3_FIND_VIRTUALENV FIRST)
set(Python3_FIND_STRATEGY LOCATION)
find_package(Python3 COMPONENTS Interpreter)
message(STATUS "Python used for the comparison: ${Python3_EXECUTABLE}")
if(Python3_FOUND AND EXISTS "${Sloth_COMPARE_SCRIPT}")
  add_test(NAME minimal_compare
           COMMAND ${Python3_EXECUTABLE} ${Sloth_COMPARE_SCRIPT}
                   -f Saves/Problem1/time_specialized.csv
                      ${CMAKE_CURRENT_SOURCE_DIR}/ref/time_specialized.csv
                   -c -1 -e relativeEpsilon -k 1e-12 -w 1)
  # Not run (instead of comparing old or partial results) if minimal_run fails
  set_tests_properties(minimal_compare PROPERTIES FIXTURES_REQUIRED minimal_results)
elseif(NOT Python3_FOUND)
  message(STATUS "Python 3 not found: comparison with the reference results disabled")
else()
  message(STATUS "Comparison script not found ('${Sloth_COMPARE_SCRIPT}'): "
                 "comparison with the reference results disabled")
endif()
