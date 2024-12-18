
create_test("AllenCahn2Dtest7" "AllenCahn2Dtest7" FALSE "2D EXE AC" 1)
create_col_comparison("CompareAllenCahn2Dtest7" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" 4 relativeEpsilon 1e-12 FALSE "AllenCahn2Dtest7" "2D AC")
create_col_comparison("CompareAllenCahn2Dtest7_diff" "time_specialized.csv" "Saves/Problem2/time_specialized.csv" 4 relativeEpsilon 1e-12 FALSE "AllenCahn2Dtest7" "2D AC")
