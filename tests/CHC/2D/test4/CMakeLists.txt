
create_test("CahnHilliard2Dtest4" "CahnHilliard2Dtest4" FALSE "2D EXE CH" 1)
create_col_comparison("CompareCahnHilliard2Dtest1" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2Dtest1" "2D AC")
