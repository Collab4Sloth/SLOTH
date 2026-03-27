create_test("HeatTransfer1Dtest1" "HeatTransfer1Dtest1" FALSE "1D EXE HEAT" 1)
create_col_comparison("CompareHeatTransfer1Dtest1" "time_specialized.csv" "Saves/Heat/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "HeatTransfer1Dtest1" "1D HEAT")
