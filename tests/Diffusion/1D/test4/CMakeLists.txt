create_test("Diffusion1Dtest4" "Diffusion1Dtest4" FALSE "1D EXE DIFFU" 4)
create_col_comparison("CompareDiffusion1Dtest4" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2Dtest2" "AC")

