
create_test("CahnHilliard3Dtest1" "CahnHilliard3Dtest1" FALSE "3D EXE CH" 1)
create_col_comparison("CompareCahnHilliard3Dtest1" "time_specialized.csv" "Saves/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "CahnHilliard3Dtest1" "3D CH")
