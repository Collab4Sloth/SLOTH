execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

create_test("HeatTransfer2Dtest8" "HeatTransfer2Dtest8" FALSE "2D EXE HEAT" 1)
create_col_comparison("CompareHeatTransfer2Dtest8" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "HeatTransfer2Dtest8" "2D HEAT")
