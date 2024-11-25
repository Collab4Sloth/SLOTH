create_test("AllenCahn2DPeriodictest3" "AllenCahn2DPeriodictest3" FALSE "2D EXE AC" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/convergence_study.py ${CMAKE_CURRENT_BINARY_DIR}/convergence_study.py COPYONLY)
create_col_comparison("CompareAllenCahn2DPeriodictest3" "time_specialized.csv" "Saves_order_1_Nx8/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "AllenCahn2DPeriodictest3" "2D AC")
create_col_comparison_convergence("CompareAllenCahn2DPeriodictest3Convergence" "convergence_output_ref.csv" "convergence_output.csv" -1 absolute 1e-16 FALSE "AllenCahn2DPeriodictest3" "2D AC" 3.141592653589793 100000)
