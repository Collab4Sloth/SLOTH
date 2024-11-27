
create_test("AllenCahn2Dtest6" "AllenCahn2Dtest6" FALSE "2D EXE AC" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/convergence_study.py ${CMAKE_CURRENT_BINARY_DIR}/convergence_study.py COPYONLY)
create_col_comparison("CompareAllenCahn2Dtest6" "time_specialized.csv" "Saves_order_2_Nx16_dt0.010000/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2Dtest6" "AC")
create_col_comparison_convergence("CompareAllenCahn2Dtest6Convergence" "convergence_output_ref.csv" "convergence_output.csv" -1 absolute 1e-16 FALSE "AllenCahn2Dtest6" "2D AC" 1. 8e-1)
