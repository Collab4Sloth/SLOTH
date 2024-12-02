create_test("Diffusion1Dtest1" "Diffusion1Dtest1" FALSE "1D EXE DIFFU" 1)
create_col_comparison("CompareDiffusion1Dtest1" "time_specialized.csv" "Saves_Nx_10/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "Diffusion1Dtest1" "1D DIFFU")



