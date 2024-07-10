
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(WARNINGS -Wall -Wextra -Wno-deprecated -Wparentheses -Wreturn-type -Wmissing-declarations -fmessage-length=0 -Wunused -Wfatal-errors -Wpointer-arith -Wcast-align -Wwrite-strings -Wctor-dtor-privacy -Wnon-virtual-dtor -Woverloaded-virtual -Wfloat-equal -Wno-endif-labels -Wsign-compare -Wmissing-format-attribute -Wno-multichar -Wno-deprecated-declarations -Wpacked -Wredundant-decls -Wdisabled-optimization -Wunknown-pragmas -Wundef -Wreorder)
set(DEBUG_OPTIONS -g -O0 -gdwarf-4 -ffpe-trap=invalid,zero,overflow,underflow -fcheck=all ${WARNINGS})

set(OPTIM_OPTIONS -g -O2 ${WARNINGS})
set(COVERAGE_OPTIONS -g -O0 --coverage)
set(COVERAGE_LINK_OPTIONS --coverage)

# COMPARISON test
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

# CREATION  test
function(create_test exe_name test_name test_will_fail test_label)
  set(CURRENT_EXE ${exe_name})

  add_executable(${CURRENT_EXE} main.cpp)
  target_link_libraries(${CURRENT_EXE} Sloth)

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

    target_compile_options(${CURRENT_EXE} PRIVATE ${COVERAGE_OPTIONS})
    link_libraries(gcov)
    set_target_properties(${CURRENT_EXE} PROPERTIES LINK_FLAGS ${COVERAGE_LINK_OPTIONS})

  ELSE()
    message("Default build : Release ")
  ENDIF()

  install(TARGETS ${CURRENT_EXE})
endfunction(create_test)
