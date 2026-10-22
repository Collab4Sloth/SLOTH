# execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
#     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

create_test("Sintering2Dtest1AMR" "Sintering2Dtest1AMR" FALSE "2D EXE CHAC AMR" 2)

create_col_comparison("CompareSintering1AMR" "time_specialized.csv" "Saves/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 5e-4 FALSE "Sintering2Dtest1AMR" "2D CHAC AMR")
