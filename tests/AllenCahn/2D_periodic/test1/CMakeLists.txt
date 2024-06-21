create_test("AllenCahn2DPeriodictest1" "AllenCahn2DPeriodictest1" FALSE "2D EXE AC")
create_col_comparison("CompareAllenCahn2DPeriodictest1_EulerExplicit" "Problem1_EulerExplicit/time_specialized.csv" "Saves/Problem1_EulerExplicit/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2DPeriodictest1" "2D AC")
create_col_comparison("CompareAllenCahn2DPeriodictest1_EulerImplicit" "Problem1_EulerImplicit/time_specialized.csv" "Saves/Problem1_EulerImplicit/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2DPeriodictest1" "2D AC")
