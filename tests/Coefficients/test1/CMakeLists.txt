set(CURRENT_EXE CoefficientTest)
add_executable(${CURRENT_EXE} main.cpp)
target_link_libraries(${CURRENT_EXE} Sloth)

add_test(NAME CoefficientTest0 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} --test-case 0)
set_tests_properties(CoefficientTest0 PROPERTIES WILL_FAIL FALSE LABELS "Coefficient")

add_test(NAME CoefficientTest1 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} --test-case 1)
set_tests_properties(CoefficientTest1 PROPERTIES WILL_FAIL TRUE LABELS "Coefficient")

add_test(NAME CoefficientTest2 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} --test-case 2)
set_tests_properties(CoefficientTest2 PROPERTIES WILL_FAIL TRUE LABELS "Coefficient")

add_test(NAME CoefficientTest3 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} --test-case 3)
set_tests_properties(CoefficientTest3 PROPERTIES WILL_FAIL FALSE LABELS "Coefficient")

add_test(NAME CoefficientTest4 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} --test-case 4)
set_tests_properties(CoefficientTest4 PROPERTIES WILL_FAIL FALSE LABELS "Coefficient")

add_test(NAME CoefficientTest5 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} --test-case 5)
set_tests_properties(CoefficientTest5 PROPERTIES WILL_FAIL FALSE LABELS "Coefficient")


set_compile_options(${CURRENT_EXE})

add_dependencies(tests ${CURRENT_EXE})

# TODO(CCI) Uncomment?
# install(TARGETS ${CURRENT_EXE} DESTINATION tests)


