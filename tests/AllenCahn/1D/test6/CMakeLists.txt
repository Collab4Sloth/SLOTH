create_test("AllenCahn1Dtest6" "AllenCahn1Dtest6" FALSE "1D EXE AC ADAPTIVE_DT" 2)
create_col_comparison("CompareAllenCahn1Dtest6" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn1Dtest6" "1D AC ADAPTIVE_DT")
