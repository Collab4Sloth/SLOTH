execute_process(COMMAND python3 ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
create_test("CahnHilliard2Dtest3" "CahnHilliard2Dtest3" FALSE "2D EXE CH" 1)
create_col_comparison("CompareCahnHilliard2Dtest3" "time_specialized.csv" "Saves/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "CahnHilliard2Dtest3" "2D CH")
