create_test("Sintering2Dtest1AMR" "Sintering2Dtest1AMR" FALSE "2D EXE CHAC AMR" 2)

create_col_comparison("CompareSintering1AMR" "time_specialized.csv" "Saves/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "Sintering2Dtest1AMR" "2D CHAC AMR")
