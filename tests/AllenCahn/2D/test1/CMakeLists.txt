
create_test("AllenCahn2Dtest1" "AllenCahn2Dtest1" FALSE "2D EXE AC")
create_col_comparison("CompareAllenCahn2Dtest1" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2Dtest1" "2D AC")
