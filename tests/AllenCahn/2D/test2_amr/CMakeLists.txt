create_test("AllenCahn2Dtest2AMR" "AllenCahn2Dtest2AMR" FALSE "EXE AC AMR" 2)

configure_file(${CMAKE_SOURCE_DIR}/Mesh-examples/camembert2D.msh ${CMAKE_CURRENT_BINARY_DIR}/camembert2D.msh COPYONLY)

create_col_comparison("CompareAllenCahn2Dtest2AMR" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2Dtest2AMR" "AC AMR")

