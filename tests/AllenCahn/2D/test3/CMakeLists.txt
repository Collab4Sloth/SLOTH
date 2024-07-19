
create_test("AllenCahn2Dtest3" "AllenCahn2Dtest3" FALSE "2D EXE AC")
create_col_comparison("CompareAllenCahn2Dtest3" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2Dtest3" "2D AC")
