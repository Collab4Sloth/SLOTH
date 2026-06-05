execute_process(COMMAND python3 ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
create_test("Grains3Dtest2" "Grains3Dtest2" FALSE "3D EXE AC" 2)
configure_file(${CMAKE_SOURCE_DIR}/tests/Studies/multigrains/test1/initial-grain.txt ${CMAKE_CURRENT_BINARY_DIR}/initial-grain.txt COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plot.py ${CMAKE_CURRENT_BINARY_DIR}/plot.py COPYONLY)

create_col_comparison("CompareGrains3D" "time_specialized.csv" "Saves/GrainsProblem/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "Grains3Dtest2" "3D AC")


add_test(
    NAME Grains3Dtest2_plot
    COMMAND python3 ${CMAKE_CURRENT_BINARY_DIR}/plot.py
)

set_tests_properties(
    Grains3Dtest2_plot
    PROPERTIES DEPENDS Grains3Dtest2
)
