
create_test("AllenCahn2Dtest1_splitting" "AllenCahn2Dtest1_splitting" FALSE "2D EXE AC" 1)
create_col_comparison("CompareAllenCahn2Dtest1_splitting" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2Dtest1_splitting" "2D AC")
