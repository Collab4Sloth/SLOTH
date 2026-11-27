execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
create_test("AllenCahn2Dtest1Variational" "AllenCahn2Dtest1Variational" FALSE "2D EXE AC" 2)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/convergence_study.py ${CMAKE_CURRENT_BINARY_DIR}/convergence_study.py COPYONLY)
create_col_comparison_convergence("CompareAllenCahn2Dtest1VariationalConvergence" "convergence_output_ref.csv" "convergence_output.csv" -1 absolute 1e-16 FALSE "AllenCahn2Dtest1Variational" "2D AC" 1. 100000)
