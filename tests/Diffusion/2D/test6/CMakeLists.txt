execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

create_test("Diffusion2Dtest6" "Diffusion2Dtest6" FALSE "2D EXE DIFFU" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/time_convergence.py ${CMAKE_CURRENT_BINARY_DIR}/time_convergence.py COPYONLY)
# Observed orders in time: EulerImplicit 1, BDF2 2, CrankNicolson (theta = 1/2) 2, with or
# without implicit Euler first step
add_test(NAME CheckDiffusion2Dtest6TimeConvergence
         COMMAND ${Python3_EXECUTABLE} time_convergence.py --check --expected "EulerImplicit:1,BDF2:2,CrankNicolson:2,CrankNicolsonStart:2"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(CheckDiffusion2Dtest6TimeConvergence PROPERTIES
                     DEPENDS Diffusion2Dtest6 LABELS "2D DIFFU")

# Non-regression on the error values
create_col_comparison("CompareDiffusion2Dtest6TimeConvergence" "convergence_output_ref.csv"
                      "time_convergence_output.csv" -1 relativeEpsilon 1e-6 FALSE
                      "CheckDiffusion2Dtest6TimeConvergence" "2D DIFFU")
