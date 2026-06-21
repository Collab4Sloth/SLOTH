
create_test("CahnHilliard2Dtest1_dt" "CahnHilliard2Dtest1_dt" FALSE "2D EXE CH" 1)
create_col_comparison("CompareCahnHilliard2Dtest1_dt" "time_specialized.csv" "Saves/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "CahnHilliard2Dtest1_dt" "2D CH")
