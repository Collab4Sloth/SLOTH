create_test("AllenCahn3Dtest1" "AllenCahn3Dtest1" FALSE "3D EXE AC" 4)
create_col_comparison("CompareAllenCahn3Dtest1a" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" 2 relativeEpsilon 1e-12 FALSE "AllenCahn3Dtest1" "3D AC")
create_col_comparison("CompareAllenCahn3Dtest1b" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" 3 relativeEpsilon 5e-3 FALSE "AllenCahn3Dtest1" "3D AC")
create_col_comparison("CompareAllenCahn3Dtest1c" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" 4 relativeEpsilon 5e-3 FALSE "AllenCahn3Dtest1" "3D AC")
