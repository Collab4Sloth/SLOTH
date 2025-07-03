
create_test("CahnHilliard2DtestLogo" "CahnHilliard2DtestLogo" FALSE "2D EXE CH" 1)
configure_file(${CMAKE_SOURCE_DIR}/kernel/Mesh-examples/slothLogo.msh ${CMAKE_CURRENT_BINARY_DIR}/slothLogo.msh COPYONLY)
create_col_comparison("CompareCahnHilliard2DtestLogo" "time_specialized.csv" "Saves/CahnHilliard/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "CahnHilliard2DtestLogo" "2D CH")
