
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(WARNINGS " -Wall  -Wextra -Wno-deprecated -Wparentheses  -Wreturn-type  -Wmissing-declarations  -fmessage-length=0  -Wunused  -Wfatal-errors  -Wpointer-arith  -Wcast-align  -Wwrite-strings  -Wctor-dtor-privacy  -Wnon-virtual-dtor  -Woverloaded-virtual -Wfloat-equal  -Wno-endif-labels  -Wsign-compare  -Wmissing-format-attribute  -Wno-multichar  -Wno-deprecated-declarations  -Wpacked  -Wredundant-decls  -Wdisabled-optimization  -Wunknown-pragmas  -Wundef  -Wreorder")
set(DEBUG_OPTIONS "-g -O0 -gdwarf-4  -ffpe-trap=invalid,zero,overflow,underflow -fcheck=all ${WARNINGS}")

# set(OPTIM_OPTIONS "-g -O2 ${WARNINGS}")
set(OPTIM_OPTIONS "-g -O2 ")
set(COVERAGE_OPTIONS " --coverage ")

# # TESTS
# file(GLOB_RECURSE CLANG_TIDY_BINARIES /usr/bin/clang-tidy*)
# FIND_PROGRAM(CLANGTIDY ${CLANG_TIDY_BINARIES})
# set(CMAKE_CXX_CLANG_TIDY "${CLANGTIDY};-checks=*,map.h,tuple.h" FORCE)

# set(CMAKE_CXX_CPPLINT "cpplint;--linelength=100 --counting=total --filter=-runtime/references,-build/header_guard,-runtime/string --recursive ${CMAKE_CURRENT_SOURCE_DIR}")
function(build_cov_target _testrunner)
  FIND_PROGRAM(GCOV_PATH gcov)
  FIND_PROGRAM(LCOV_PATH lcov)
  FIND_PROGRAM(GENHTML_PATH genhtml)

  IF(NOT GCOV_PATH)
    MESSAGE(FATAL_ERROR "gcov not found! Aborting...")
  ENDIF() # NOT GCOV_PATH

  IF(NOT LCOV_PATH)
    MESSAGE(FATAL_ERROR "lcov not found! Aborting...")
  ENDIF() # NOT LCOV_PATH

  IF(NOT GENHTML_PATH)
    MESSAGE(FATAL_ERROR "genhtml not found! Aborting...")
  ENDIF() # NOT GENHTML_PATH

  ADD_CUSTOM_TARGET(${_testrunner}_coverage

    # Cleanup lcov
    ${LCOV_PATH} --directory . --zerocounters
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/${_testrunner}.info ${CMAKE_BINARY_DIR}/${_testrunner}.info.cleaned

    # Run tests
    COMMAND ctest -R ${_testrunner}

    # Capturing lcov counters and generating report
    COMMAND ${LCOV_PATH} --directory . --capture --output-file ${CMAKE_BINARY_DIR}/${_testrunner}.info
    COMMAND ${LCOV_PATH} --directory . --remove ${CMAKE_BINARY_DIR}/${_testrunner}.info '*/mfem/*' '*/Tests/*' '/usr/*' --output-file ${CMAKE_BINARY_DIR}/${_testrunner}.info.cleaned
    COMMAND ${GENHTML_PATH} -o ${CMAKE_BINARY_DIR}/Coverage_${_testrunner} ${CMAKE_BINARY_DIR}/${_testrunner}.info.cleaned

    COMMAND if [ -f '${CMAKE_BINARY_DIR}/merged.info' ] \; then ${LCOV_PATH} --add-tracefile ${CMAKE_BINARY_DIR}/${_testrunner}.info.cleaned --add-tracefile ${CMAKE_BINARY_DIR}/merged.info -o ${CMAKE_BINARY_DIR}/merged.info \; fi \;
    COMMAND if [ -f '${CMAKE_BINARY_DIR}/merged.info' ] \; then ${GENHTML_PATH} -o ${CMAKE_BINARY_DIR}/Coverage ${CMAKE_BINARY_DIR}/merged.info \; fi \;
    COMMAND if [ ! -f '${CMAKE_BINARY_DIR}/merged.info' ] \; then cp ${CMAKE_BINARY_DIR}/${_testrunner}.info.cleaned ${CMAKE_BINARY_DIR}/merged.info \; fi \;

    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Resetting code coverage counters to zero.\nProcessing code coverage counters and generating report."
  )
endfunction() # SETUP_TARGET_FOR_COVERAGE

function(create_col_comparison test_name reference_file results_file cols criterion threshold test_will_fail test_depend test_label)
  set(DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ref_${test_name})
  file(GLOB_RECURSE REF_FILE ${CMAKE_CURRENT_SOURCE_DIR}/ref/${reference_file})
  file(MAKE_DIRECTORY ${DEST_DIR})
  configure_file(${REF_FILE} ${DEST_DIR}/${reference_file} COPYONLY)
  configure_file(${REF_FILE} ${DEST_DIR}/${reference_file} COPYONLY)
  configure_file(${CMAKE_SOURCE_DIR}/tests/tools/col_compare.py ${CMAKE_CURRENT_BINARY_DIR}/col_compare.py COPYONLY)

  add_test(NAME ${test_name} COMMAND python3 col_compare.py -f ${results_file} ${DEST_DIR}/${reference_file} -c ${cols} -e ${criterion} -k ${threshold} -w 1)
  set_tests_properties(PROPERTIES WILL_FAIL ${test_will_fail})
  set_tests_properties(${test_name} PROPERTIES DEPENDS ${test_depend})
  set_tests_properties(${test_name} PROPERTIES LABELS ${test_label})
endfunction() # COMPARISON test

function(create_test exe_name test_name test_will_fail test_label)
  set(CURRENT_EXE ${exe_name})

  add_executable(${CURRENT_EXE} main.cpp)
  target_link_libraries(${CURRENT_EXE} Sloth)

  # target_include_directories(${CURRENT_EXE} PUBLIC ${MFEM_DIR}/include)
  # target_link_directories(${CURRENT_EXE} PUBLIC ${MFEM_DIR}/lib $ENV{HYPRE_DIR}/lib $ENV{MPI_DIR}/lib $ENV{PETSC_DIR}/lib)
  # target_link_libraries(${CURRENT_EXE} mfem mpi z HYPRE SuiteSparse::suitesparse metis petsc $ENV{FILESYSTEM_VAR} Sloth)
  add_test(NAME ${test_name} COMMAND ${CURRENT_EXE})
  set_tests_properties(PROPERTIES WILL_FAIL ${test_will_fail})
  set_tests_properties(${test_name} PROPERTIES LABELS ${test_label})

  IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug build.")
    set_target_properties(${CURRENT_EXE} PROPERTIES COMPILE_FLAGS ${DEBUG_OPTIONS})

  ELSEIF(CMAKE_BUILD_TYPE MATCHES Optim)
    message("Optim build.")
    set_target_properties(${CURRENT_EXE} PROPERTIES COMPILE_FLAGS ${OPTIM_OPTIONS})

  ELSEIF(CMAKE_BUILD_TYPE MATCHES Coverage)
    message("Coverage build.")
    build_cov_target(${CURRENT_EXE})

    target_link_libraries(${CURRENT_EXE} gcov)
    set_target_properties(${CURRENT_EXE} PROPERTIES COMPILE_FLAGS ${COVERAGE_OPTIONS})
    set_target_properties(${CURRENT_EXE} PROPERTIES LINK_FLAGS ${COVERAGE_OPTIONS})

  ELSE()
    message("Default build : Release ")
  ENDIF()

  install(TARGETS ${CURRENT_EXE})
endfunction(create_test)
