execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

create_test("HeatTransfer2Dtest10" "HeatTransfer2Dtest10" FALSE "2D EXE HEAT" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/time_convergence.py ${CMAKE_CURRENT_BINARY_DIR}/time_convergence.py COPYONLY)
# Observed orders in time: EulerImplicit 1, BDF2 2, CrankNicolson (theta = 1/2) 2
add_test(NAME CheckHeatTransfer2Dtest10TimeConvergence
         COMMAND ${Python3_EXECUTABLE} time_convergence.py --check --expected "EulerImplicit:1,BDF2:2,CrankNicolson:2"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(CheckHeatTransfer2Dtest10TimeConvergence PROPERTIES
                     DEPENDS HeatTransfer2Dtest10 LABELS "2D HEAT")

# Non-regression on the error values
create_col_comparison("CompareHeatTransfer2Dtest10TimeConvergence" "convergence_output_ref.csv"
                      "time_convergence_output.csv" -1 relativeEpsilon 1e-6 FALSE
                      "CheckHeatTransfer2Dtest10TimeConvergence" "2D HEAT")
