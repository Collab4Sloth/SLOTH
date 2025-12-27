
create_test("CahnHilliard2Dtest4" "CahnHilliard2Dtest4" FALSE "2D EXE CH" 2)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/convergence_study.py ${CMAKE_CURRENT_BINARY_DIR}/convergence_study.py COPYONLY)
create_col_comparison_convergence("CompareCahnHilliard2Dtest4Convergence" "convergence_output_ref.csv" "convergence_output.csv" -1 absolute 1e-16 FALSE "CahnHilliard2Dtest4" "2D CH" 1. 100000)
