execute_process(COMMAND python3 ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set(CURRENT_EXE Ostwald2Dtest1)
add_executable(${CURRENT_EXE} main.cpp)
target_link_libraries(${CURRENT_EXE} Sloth)

add_test(
    NAME Ostwald2Dtest1_base
    COMMAND mpirun -np 2 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} -t 0.05
)

add_test(
    NAME Ostwald2Dtest1_restart
    COMMAND mpirun -np 2 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} -r -s 5 -i 0.05
)

set_tests_properties(
    Ostwald2Dtest1_restart
    PROPERTIES DEPENDS Ostwald2Dtest1_base
)



create_col_comparison("CompareSaveRestart" "time_specialized.csv" "Saves/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "Ostwald2Dtest1_restart" "2D AC")
