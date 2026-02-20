execute_process(COMMAND python3 ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/CoefficientExample.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set(CURRENT_EXE CoefficientTest2)
add_executable(${CURRENT_EXE} main.cpp)
target_link_libraries(${CURRENT_EXE} Sloth)
target_include_directories(${CURRENT_EXE} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(${CURRENT_EXE} PRIVATE ${CMAKE_SOURCE_DIR}/kernel/Coefficients)


add_test(NAME CoefficientTest20 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} --test-case 0)
set_tests_properties(CoefficientTest20 PROPERTIES WILL_FAIL FALSE LABELS "Coefficient")

add_test(NAME CoefficientTest21 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} --test-case 1)
set_tests_properties(CoefficientTest21 PROPERTIES WILL_FAIL FALSE LABELS "Coefficient")


add_test(NAME CoefficientTest22 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} --test-case 2)
set_tests_properties(CoefficientTest22 PROPERTIES WILL_FAIL FALSE LABELS "Coefficient")

add_test(NAME CoefficientTest23 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} --test-case 3)
set_tests_properties(CoefficientTest23 PROPERTIES WILL_FAIL FALSE LABELS "Coefficient")

add_test(NAME CoefficientTest24 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} --test-case 4)
set_tests_properties(CoefficientTest24 PROPERTIES WILL_FAIL FALSE LABELS "Coefficient")

add_test(NAME CoefficientTest25 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} --test-case 5)
set_tests_properties(CoefficientTest25 PROPERTIES WILL_FAIL FALSE LABELS "Coefficient")

set_compile_options(${CURRENT_EXE})

