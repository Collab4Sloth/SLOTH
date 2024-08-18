create_test("Diffusion2Dtest1" "Diffusion2Dtest1" FALSE "2D EXE DIFFU" 1)

configure_file(${CMAKE_SOURCE_DIR}/kernel/Mesh-examples/star2D.msh ${CMAKE_CURRENT_BINARY_DIR}/star2D.msh COPYONLY)
create_col_comparison("CompareDiffusion2Dtest1" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "Diffusion2Dtest1" "2D DIFFU")
