create_test("HeatTransfer1Dtest3" "HeatTransfer1Dtest3" FALSE "2D EXE HEAT" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/convergence_study.py ${CMAKE_CURRENT_BINARY_DIR}/convergence_study.py COPYONLY)
create_col_comparison_convergence("CompareHeatTransfer1Dtest3Convergence" "convergence_output_ref.csv" "convergence_output.csv" -1 absolute 1e-16 FALSE "HeatTransfer1Dtest3" "1D HEAT" 0.00465 100000)



