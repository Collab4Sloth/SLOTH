
create_test("CahnHilliard2Dtest2" "CahnHilliard2Dtest2" FALSE "2D EXE CH" 1)
create_col_comparison("CompareCahnHilliard2Dtest2" "time_specialized.csv" "Saves/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "CahnHilliard2Dtest2" "2D CH")
