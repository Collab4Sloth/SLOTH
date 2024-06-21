create_test("AllenCahn3Dtest2" "AllenCahn3Dtest2" FALSE "3D EXE AC")

configure_file(${CMAKE_SOURCE_DIR}/kernel/Mesh-examples/camembert3D.msh ${CMAKE_CURRENT_BINARY_DIR}/camembert3D.msh COPYONLY)
create_col_comparison("CompareAllenCahn3Dtest2" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn3Dtest2" "3D AC")

