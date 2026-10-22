# ##############################################################################
# SLOTH unit tests (Catch2 v3)
# ##############################################################################
# Catch2 v3 is an optional dependency: without it, SLOTH and the other tests
# are built as usual and the unit tests are simply skipped.
#   macOS : brew install catch2
#   spack : spack install catch2   (then spack load catch2)
#   other : -DCatch2_DIR=<prefix>/lib/cmake/Catch2
#
# Usage:
#   make unit_tests                 # build all unit tests
#   ctest -L Unit                   # run all unit tests
#   ctest -L Parameters             # run one module
# ##############################################################################
find_package(Catch2 3 QUIET)
if(NOT Catch2_FOUND)
  message(STATUS "Catch2 v3 not found: SLOTH unit tests are disabled")
  return()
endif()
message(STATUS "Found Catch2 ${Catch2_VERSION}: SLOTH unit tests are enabled")

# ------------------------------------------------------------------------------
# Common main (MPI + Hypre initialisation, result reduced over the ranks).
# OBJECT library: compiled once, its object file is added to every executable
# that links it. Catch2::Catch2 (and not Catch2::Catch2WithMain): the main is ours.
# ------------------------------------------------------------------------------
add_library(sloth_unit_main OBJECT framework/sloth_unit_main.cpp)
target_link_libraries(sloth_unit_main PUBLIC Sloth Catch2::Catch2)
# Options of the build type (Release, Debug, Coverage...), as for the other tests
set_compile_options(sloth_unit_main)
# Tests include "sloth_unit.hpp" (helpers)
target_include_directories(sloth_unit_main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/framework)

# Target building all unit tests (and only them)
add_custom_target(unit_tests COMMENT "Build SLOTH unit tests")

# ------------------------------------------------------------------------------
# sloth_add_unit_test(NAME    <module>
#                     SOURCES <file.cpp> [<file.cpp> ...]
#                     [NPROCS <n> [<n> ...]]       # default: 1 2 (CI: at most 2)
#                     [COEFFICIENTS_JSON <file.json>])
#
# COEFFICIENTS_JSON: coefficients generated at configure time by
# scripts/GenerateCoefficient.py (as for the other SLOTH tests); the generated
# headers are written in the build directory and can be included by the tests.
#
# Creates:
#   - the executable Unit_<module> (common main + SOURCES);
#   - one ctest entry per number of ranks: Unit.<module>.np<n>,
#     run with mpirun, labels "Unit;<module>".
# ------------------------------------------------------------------------------
function(sloth_add_unit_test)
  # Reads the keyword arguments into UT_NAME, UT_SOURCES, UT_NPROCS
  cmake_parse_arguments(UT "" "NAME;COEFFICIENTS_JSON" "SOURCES;NPROCS" ${ARGN})
  if(NOT UT_NAME OR NOT UT_SOURCES)
    message(FATAL_ERROR "sloth_add_unit_test: NAME and SOURCES are required")
  endif()
  if(NOT UT_NPROCS)
    set(UT_NPROCS 1 2)
  endif()

  # Executable
  set(exe Unit_${UT_NAME})
  add_executable(${exe} ${UT_SOURCES})
  target_link_libraries(${exe} PRIVATE sloth_unit_main)
  # Options of the build type (Release, Debug, Coverage...), as for the other
  # tests. In Coverage, this instruments the kernel/include code compiled in
  # the test (templates, inline functions).
  set_compile_options(${exe})

  # Coefficients generated from JSON (optional)
  if(UT_COEFFICIENTS_JSON)
    set(gen_dir ${CMAKE_CURRENT_BINARY_DIR}/generated_${UT_NAME})
    file(MAKE_DIRECTORY ${gen_dir})
    execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py
                            -r -f ${CMAKE_CURRENT_SOURCE_DIR}/${UT_COEFFICIENTS_JSON}
                    WORKING_DIRECTORY ${gen_dir}
                    RESULT_VARIABLE gen_result)
    if(NOT gen_result EQUAL 0)
      message(FATAL_ERROR "sloth_add_unit_test(${UT_NAME}): coefficient generation failed")
    endif()
    target_include_directories(${exe} PRIVATE ${gen_dir})
    # Regenerate when the JSON file changes
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                 ${CMAKE_CURRENT_SOURCE_DIR}/${UT_COEFFICIENTS_JSON})
  endif()

  # ctest entries
  foreach(np IN LISTS UT_NPROCS)
    set(test_name Unit.${UT_NAME}.np${np})
    add_test(NAME ${test_name}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS}
                     $<TARGET_FILE:${exe}> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(${test_name} PROPERTIES
                         LABELS "Unit;${UT_NAME}"
                         PROCESSORS ${np}
                         # A unit test runs in seconds: a hang (MPI deadlock...)
                         # fails after 2 minutes instead of blocking the CI
                         TIMEOUT 120)
  endforeach()

  # Built by "make unit_tests" and by "make tests"
  add_dependencies(unit_tests ${exe})
  add_dependencies(tests ${exe})
endfunction()

# ##############################################################################
# Modules
# ##############################################################################
sloth_add_unit_test(NAME AMR
                    SOURCES AMR/test_SlothErrorEstimators.cpp)

sloth_add_unit_test(NAME AnalyticalFunctions
                    SOURCES AnalyticalFunctions/test_AnalyticalFunctions.cpp)

sloth_add_unit_test(NAME BCs
                    SOURCES BCs/test_Boundary.cpp
                    NPROCS 1)

sloth_add_unit_test(NAME Coefficients
                    SOURCES Coefficients/test_FunctionCoefficients.cpp
                            Coefficients/test_SlothBaseCoefficient.cpp
                            Coefficients/test_Coefficients.cpp
                            Coefficients/test_GeneratedCoefficients.cpp
                            Coefficients/test_FileCoefficient.cpp
                            Coefficients/test_MfemCoefficients.cpp
                    COEFFICIENTS_JSON Coefficients/CoefficientExample.json)

sloth_add_unit_test(NAME Convergence
                    SOURCES Convergence/test_Convergence.cpp
                    NPROCS 1)

sloth_add_unit_test(NAME Framework
                    SOURCES framework/test_mpi.cpp)

sloth_add_unit_test(NAME Glossary
                    SOURCES Glossary/test_Glossary.cpp
                    NPROCS 1)

sloth_add_unit_test(NAME Integrators
                    SOURCES Integrators/test_SlothGridFunction.cpp)

sloth_add_unit_test(NAME Operators
                    SOURCES Operators/test_ReducedOperators.cpp)

sloth_add_unit_test(NAME Options
                    SOURCES Options/test_Options.cpp
                    NPROCS 1)

sloth_add_unit_test(NAME Parameters
                    SOURCES Parameters/test_Parameter.cpp
                            Parameters/test_Parameters.cpp)

sloth_add_unit_test(NAME PreProcessing
                    SOURCES PreProcessing/test_TabulatedFunction.cpp
                    NPROCS 1)

sloth_add_unit_test(NAME Problems
                    SOURCES Problems/test_CahnHilliardDirichlet.cpp)

sloth_add_unit_test(NAME Property
                    SOURCES Property/test_InterDiffusionCoefficient.cpp
                    NPROCS 1)

sloth_add_unit_test(NAME Solvers
                    SOURCES Solvers/test_Solvers.cpp)

sloth_add_unit_test(NAME Utils
                    SOURCES Utils/test_Utils.cpp
                    NPROCS 1)

# ##############################################################################
# Coverage of the unit tests alone (CMAKE_BUILD_TYPE=Coverage)
# ##############################################################################
#   make unit_coverage
# 1. resets the coverage counters (*.gcda) of the whole build directory;
# 2. runs the unit tests only (ctest -L Unit);
# 3. gcovr: searches the whole build directory, so that the code of
#    kernel/include compiled in the test executables (templates, inline
#    functions) is taken into account, and writes
#      Coverage/unit/coverage.html  (detailed report, file by file)
#      Coverage/unit/coverage.json  (input of the script below);
# 4. tools/uncovered_functions.py: number of uncovered functions per module.
#    Details for one module:
#      python3 tests/Unit/tools/uncovered_functions.py \
#              build/Coverage/unit/coverage.json -m Parameters
# ##############################################################################
if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
  set(UNIT_COV_DIR ${CMAKE_BINARY_DIR}/Coverage/unit)
  # With Clang (macOS), the gcov-compatible tool is llvm-cov gcov
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(APPLE)
      set(UNIT_GCOV_OPTION --gcov-executable "xcrun llvm-cov gcov")
    else()
      set(UNIT_GCOV_OPTION --gcov-executable "llvm-cov gcov")
    endif()
  endif()

  add_custom_target(unit_coverage
    COMMAND find ${CMAKE_BINARY_DIR} -name "*.gcda" -delete
    COMMAND ${CMAKE_CTEST_COMMAND} -L Unit
    COMMAND ${CMAKE_COMMAND} -E make_directory ${UNIT_COV_DIR}
    COMMAND gcovr -r ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} ${UNIT_GCOV_OPTION}
            --filter ${CMAKE_SOURCE_DIR}/kernel/src
            --filter ${CMAKE_SOURCE_DIR}/kernel/include
            --exclude "${CMAKE_SOURCE_DIR}/kernel/include/MAToolsProfiling/.*"
            --exclude-unreachable-branches
            --json ${UNIT_COV_DIR}/coverage.json
            --html-details ${UNIT_COV_DIR}/coverage.html
            --print-summary
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/uncovered_functions.py
            ${UNIT_COV_DIR}/coverage.json
    DEPENDS unit_tests
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Coverage of the unit tests (see Coverage/unit/coverage.html)"
    VERBATIM)
endif()
