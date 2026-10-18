execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

create_test("CahnHilliard2Dtest6" "CahnHilliard2Dtest6" FALSE "2D EXE CH ADVECTION AMR" 2)
create_col_comparison("CompareCahnHilliard2Dtest6_periodic" "periodic_time_specialized.csv" "Saves_periodic/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "CahnHilliard2Dtest6_periodic" "2D CH")
create_col_comparison("CompareCahnHilliard2Dtest6_open" "open_time_specialized.csv" "Saves_open/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "CahnHilliard2Dtest6_open" "2D CH")


set_tests_properties(
    CompareCahnHilliard2Dtest6_open
    PROPERTIES DEPENDS CahnHilliard2Dtest6
)
set_tests_properties(
    CompareCahnHilliard2Dtest6_periodic
    PROPERTIES DEPENDS CahnHilliard2Dtest6
)