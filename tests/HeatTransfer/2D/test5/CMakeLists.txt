execute_process(COMMAND python3 ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

create_test("HeatTransfer2Dtest5" "HeatTransfer2Dtest5" FALSE "2D EXE HEAT" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/convergence_study.py ${CMAKE_CURRENT_BINARY_DIR}/convergence_study.py COPYONLY)
create_col_comparison_convergence("CompareHeat2Dtest5Convergence" "convergence_output_ref.csv" "convergence_output.csv" -1 absolute 1e-16 FALSE "HeatTransfer2Dtest5" "2D Heat" 1. 100000)