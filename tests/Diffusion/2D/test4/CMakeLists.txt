create_test("Diffusion2Dtest4" "Diffusion2Dtest4" FALSE "2D EXE DIFFU" 2)

create_col_comparison("CompareDiffusion2Dtest4_inter_stab" "InterDiffusion/time_specialized.csv" "Saves/InterDiffusion/time_specialized.csv" -1 relativeEpsilon 1e-4 FALSE "Diffusion2Dtest4" "2D DIFFU")
create_col_comparison("CompareDiffusion2Dtest4_fick" "Diffusion/time_specialized.csv" "Saves/Diffusion/time_specialized.csv" -1 relativeEpsilon 1e-4 FALSE "Diffusion2Dtest4" "2D DIFFU")


