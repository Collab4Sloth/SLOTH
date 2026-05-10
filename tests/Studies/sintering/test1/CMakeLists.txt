# execute_process(COMMAND python3 ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
#     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
create_test("Sintering2Dtest1" "Sintering2Dtest1" FALSE "2D EXE CHAC" 2)

create_col_comparison("CompareSintering1" "time_specialized.csv" "Saves/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "Sintering2Dtest1" "2D CHAC")
