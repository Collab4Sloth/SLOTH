set(CURRENT_EXE CoefficientTest3)
add_executable(${CURRENT_EXE} main.cpp)
target_link_libraries(${CURRENT_EXE} Sloth)
target_include_directories(${CURRENT_EXE} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME CoefficientTest3 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} --test-case 0)
set_tests_properties(CoefficientTest3 PROPERTIES WILL_FAIL FALSE LABELS "Coefficient")

set_compile_options(${CURRENT_EXE})
add_dependencies(tests ${CURRENT_EXE})
