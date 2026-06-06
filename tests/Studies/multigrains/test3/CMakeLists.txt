execute_process(COMMAND python3 ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
create_test("Grains2Dtest3" "Grains2Dtest3" FALSE "2D EXE AC" 2)
configure_file(${CMAKE_SOURCE_DIR}/tests/Studies/multigrains/test3/initial-grain.txt ${CMAKE_CURRENT_BINARY_DIR}/initial-grain.txt COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plot.py ${CMAKE_CURRENT_BINARY_DIR}/plot.py COPYONLY)

create_col_comparison("CompareGrains" "time_specialized.csv" "Saves/GrainsProblem/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "Grains2Dtest3" "2D AC")


add_test(
    NAME Grains2Dtest3_plot
    COMMAND python3 ${CMAKE_CURRENT_BINARY_DIR}/plot.py
)

set_tests_properties(
    Grains2Dtest3_plot
    PROPERTIES DEPENDS Grains2Dtest3
)
