
create_test("AllenCahn2DPeriodictest2" "AllenCahn2DPeriodictest2" FALSE "2D EXE AC")
create_col_comparison("CompareAllenCahn2DPeriodictest2" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2DPeriodictest2" "2D AC")



