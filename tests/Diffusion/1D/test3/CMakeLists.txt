create_test("Diffusion1Dtest3" "Diffusion1Dtest3" FALSE "1D EXE DIFFU" 4)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/convergence_study.py ${CMAKE_CURRENT_BINARY_DIR}/convergence_study.py COPYONLY)
create_col_comparison_convergence("CompareDiffusion1Dtest3Convergence" "convergence_output_ref.csv" "convergence_output.csv" -1 absolute 1e-16 FALSE "Diffusion1Dtest3" "1D EXE DIFFU" 1.e-3 100000)



