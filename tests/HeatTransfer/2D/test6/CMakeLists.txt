execute_process(COMMAND python3 ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/dirichlet.csv ${CMAKE_CURRENT_BINARY_DIR}/dirichlet.csv COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/neumann.csv ${CMAKE_CURRENT_BINARY_DIR}/neumann.csv COPYONLY)
create_test("HeatTransfer2Dtest6" "HeatTransfer2Dtest6" FALSE "2D EXE HEAT" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/convergence_study.py ${CMAKE_CURRENT_BINARY_DIR}/convergence_study.py COPYONLY)
create_col_comparison_convergence("CompareHeat2Dtest6Convergence" "convergence_output_ref.csv" "convergence_output.csv" -1 absolute 1e-16 FALSE "HeatTransfer2Dtest6" "2D Heat" 1. 100000)