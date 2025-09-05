
create_test("AllenCahn2Dtest1" "AllenCahn2Dtest1" FALSE "2D EXE AC" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/convergence_study.py ${CMAKE_CURRENT_BINARY_DIR}/convergence_study.py COPYONLY)
create_col_comparison_convergence("CompareAllenCahn2Dtest1Convergence" "convergence_output_ref.csv" "convergence_output.csv" -1 absolute 1e-16 FALSE "AllenCahn2Dtest1" "2D AC" 1. 100000)
