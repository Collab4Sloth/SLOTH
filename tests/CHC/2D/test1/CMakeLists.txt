
create_test("CahnHilliard2Dtest1" "CahnHilliard2Dtest1" FALSE "2D EXE AC" 1)
create_col_comparison("CompareCahnHilliard2Dtest1" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2Dtest1" "2D AC")
