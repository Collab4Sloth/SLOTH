execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

create_test("CahnHilliard3Dtest1" "CahnHilliard3Dtest1" FALSE "3D EXE CH ADVECTION AMR" 2)



create_col_comparison("CompareCahnHilliard3Dtest1_periodic" "periodic_time_specialized.csv" "Saves_periodic/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "CahnHilliard3Dtest1_periodic" "3D CH")
create_col_comparison("CompareCahnHilliard3Dtest1_open" "open_time_specialized.csv" "Saves_open/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "CahnHilliard3Dtest1_open" "3D CH")


set_tests_properties(
    CompareCahnHilliard3Dtest1_open
    PROPERTIES DEPENDS CahnHilliard3Dtest1
)
set_tests_properties(
    CompareCahnHilliard3Dtest1_periodic
    PROPERTIES DEPENDS CahnHilliard3Dtest1
)