execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


create_test("AllenCahn2Dtest3" "AllenCahn2Dtest3" FALSE "2D EXE AC" 1)
create_col_comparison("CompareAllenCahn2Dtest3" "time_specialized.csv" "Saves_order_1_Nx4/Problem1/time_specialized.csv" -1 absolute 1e-15 FALSE "AllenCahn2Dtest3" "2D AC")
