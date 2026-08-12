execute_process(COMMAND python3 ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
create_test("CahnHilliard2Dtest2AMR" "CahnHilliard2Dtest2AMR" FALSE "2D EXE CH AMR" 1)
create_col_comparison("CompareCahnHilliard2Dtest2AMR" "time_specialized.csv" "Saves/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "CahnHilliard2Dtest2AMR" "2D CH AMR")
