create_test("Diffusion2Dtest3" "Diffusion2Dtest3" FALSE "2D EXE DIFFU" 1)
create_col_comparison("CompareDiffusion2Dtest3" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "Diffusion2Dtest3" "2D DIFFU")
