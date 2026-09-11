execute_process(COMMAND python3 ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
create_test("AMRGrains2Dtest1" "AMRGrains2Dtest1" FALSE "2D EXE AC" 2)
configure_file(${CMAKE_SOURCE_DIR}/tests/Studies/multigrains/test1_amr/initial-grain.txt ${CMAKE_CURRENT_BINARY_DIR}/initial-grain.txt COPYONLY)


create_col_comparison("CompareGrainsAMR" "time_specialized.csv" "Saves/GrainsProblem/Default_NonLinear_problem29_time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AMRGrains2Dtest1" "2D AC")
