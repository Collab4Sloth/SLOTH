create_test("AllenCahn1Dtest1" "AllenCahn1Dtest1" FALSE "1D EXE AC")
create_col_comparison("CompareAllenCahn1Dtest1Problem1" "Problem1/time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn1Dtest1" "1D AC")
create_col_comparison("CompareAllenCahn1Dtest1Problem2" "Problem2/time_specialized.csv" "Saves/Problem2/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn1Dtest1" "1D AC")
create_col_comparison("CompareAllenCahn1Dtest1Problem3" "Problem3/time_specialized.csv" "Saves/Problem3/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn1Dtest1" "1D AC")



