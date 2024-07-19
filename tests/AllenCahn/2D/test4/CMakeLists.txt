
create_test("AllenCahn2Dtest4" "AllenCahn2Dtest4" FALSE "2D EXE AC")
create_col_comparison("CompareAllenCahn2Dtest4" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2Dtest4" "2D AC")
