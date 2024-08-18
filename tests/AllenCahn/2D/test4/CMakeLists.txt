
create_test("AllenCahn2Dtest4" "AllenCahn2Dtest4" FALSE "2D EXE AC" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/convergence_study.py ${CMAKE_CURRENT_BINARY_DIR}/convergence_study.py COPYONLY)
create_col_comparison("CompareAllenCahn2Dtest4" "time_specialized.csv" "Saves_order_1_Nx8/Problem1/time_specialized.csv" -1 absolute 1e-16 FALSE "AllenCahn2Dtest4" "2D AC")
create_col_comparison_convergence("CompareAllenCahn2Dtest4Convergence" "convergence_output_ref.csv" "convergence_output.csv" -1 absolute 1e-16 FALSE "AllenCahn2Dtest4" "2D AC" 1.)
