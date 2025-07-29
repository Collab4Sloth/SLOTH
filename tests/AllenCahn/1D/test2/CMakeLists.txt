create_test("AllenCahn1Dtest2" "AllenCahn1Dtest2" FALSE "1D EXE AC" 1)
create_col_comparison("CompareAllenCahn1Dtest2" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn1Dtest2" "1D AC")
create_col_comparison("CompareAllenCahn1Dtest2Iso" "phi1_iso_computation.csv" "Saves/Problem1/phi1_iso_computation.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn1Dtest2" "1D AC")



