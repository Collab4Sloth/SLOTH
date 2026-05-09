
create_test("Nucleation1" "Nucleation1" FALSE "2D EXE AC" 2)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plot.py ${CMAKE_CURRENT_BINARY_DIR}/plot.py COPYONLY)

create_col_comparison("CompareNucleation1_099" "ref099/time_specialized.csv" "Saves_099/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "Nucleation1" "1D AC")
create_col_comparison("CompareNucleation1_1" "ref1/time_specialized.csv" "Saves_1/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "Nucleation1" "1D AC")
create_col_comparison("CompareNucleation1_101" "ref101/time_specialized.csv" "Saves_101/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "Nucleation1" "1D AC")



add_test(
    NAME Nucleation1_plot
    COMMAND python3 ${CMAKE_CURRENT_BINARY_DIR}/plot.py
)

set_tests_properties(
    Nucleation1_plot
    PROPERTIES DEPENDS Nucleation1
)

