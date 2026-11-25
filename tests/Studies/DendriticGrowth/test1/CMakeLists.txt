execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
create_test("DendriticGrowth2Dtest1" "DendriticGrowth2Dtest1" FALSE "2D EXE AC AMR" 2)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plot.py ${CMAKE_CURRENT_BINARY_DIR}/plot.py COPYONLY)
create_col_comparison("CompareDendriticGrowth2Dtest1" "time_specialized.csv" "Saves/PhaseField/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "DendriticGrowth2Dtest1" "2D AC AMR")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pfhub_bm3_references.csv ${CMAKE_CURRENT_BINARY_DIR}/pfhub_bm3_references.csv COPYONLY)
