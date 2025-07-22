
create_test("CahnHilliard2Dtest4" "CahnHilliard2Dtest4" FALSE "2D EXE CH" 1)
create_col_comparison("CompareCahnHilliard2Dtest2" "time_specialized.csv" "Saves/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "CahnHilliard2Dtest4" "2D CH")
