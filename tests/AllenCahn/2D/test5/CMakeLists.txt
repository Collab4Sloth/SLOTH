
create_test("AllenCahn2Dtest5" "AllenCahn2Dtest5" FALSE "2D EXE AC" 1)
create_col_comparison("CompareAllenCahn2Dtest5" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2Dtest5" "2D AC")
