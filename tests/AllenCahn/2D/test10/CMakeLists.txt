create_test("AllenCahn2Dtest10" "AllenCahn2Dtest10" FALSE "EXE AC" 2)

configure_file(${CMAKE_SOURCE_DIR}/kernel/Mesh-examples/pellet2Dinclusion.msh ${CMAKE_CURRENT_BINARY_DIR}/pellet2Dinclusion.msh COPYONLY)

create_col_comparison("CompareAllenCahn2Dtest10" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 absolute 2e-5 FALSE "AllenCahn2Dtest10" "AC")

