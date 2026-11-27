execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


create_test("CahnHilliard2Dtest5_bdf2" "CahnHilliard2Dtest5_bdf2" FALSE "2D EXE CH" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/time_convergence.py ${CMAKE_CURRENT_BINARY_DIR}/time_convergence.py COPYONLY)
# Observed orders in time: EulerImplicit 1, BDF2 2
add_test(NAME CheckCahnHilliard2Dtest5_bdf2TimeConvergence
         COMMAND ${Python3_EXECUTABLE} time_convergence.py --check
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(CheckCahnHilliard2Dtest5_bdf2TimeConvergence PROPERTIES
                     DEPENDS CahnHilliard2Dtest5_bdf2 LABELS "2D CH")
                     
# Non-regression on the error values 
create_col_comparison("CompareCahnHilliard2Dtest5_bdf2TimeConvergence" "convergence_output_ref.csv"
                      "time_convergence_output.csv" -1 relativeEpsilon 1e-6 FALSE
                      "CheckCahnHilliard2Dtest5_bdf2TimeConvergence" "2D CH")
