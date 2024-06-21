create_test("AllenCahn3Dtest1" "AllenCahn3Dtest1" FALSE "3D EXE AC")
create_col_comparison("CompareAllenCahn3Dtest1" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn3Dtest1" "3D AC")
