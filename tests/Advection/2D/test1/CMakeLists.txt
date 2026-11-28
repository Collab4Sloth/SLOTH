execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

create_test("Advection2Dtest1" "Advection2Dtest1" FALSE "2D EXE ADVECTION" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/time_convergence.py ${CMAKE_CURRENT_BINARY_DIR}/time_convergence.py COPYONLY)
# Observed orders in time: EulerImplicit 1, CrankNicolson (theta = 1/2) 2
add_test(NAME CheckAdvection2Dtest1TimeConvergence
         COMMAND ${Python3_EXECUTABLE} time_convergence.py --check --expected "EulerImplicit:1,CrankNicolson:2"
                 --label "exact transport"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(CheckAdvection2Dtest1TimeConvergence PROPERTIES
                     DEPENDS Advection2Dtest1 LABELS "2D ADVECTION")

# Non-regression on the error values
create_col_comparison("CompareAdvection2Dtest1TimeConvergence" "convergence_output_ref.csv"
                      "time_convergence_output.csv" -1 relativeEpsilon 1e-6 FALSE
                      "CheckAdvection2Dtest1TimeConvergence" "2D ADVECTION")
