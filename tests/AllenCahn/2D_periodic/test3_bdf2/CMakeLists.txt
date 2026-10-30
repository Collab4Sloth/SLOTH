execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

create_test("AllenCahn2DPeriodictest3_bdf2" "AllenCahn2DPeriodictest3_bdf2" FALSE "2D EXE AC" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/time_convergence.py ${CMAKE_CURRENT_BINARY_DIR}/time_convergence.py COPYONLY)
# Observed orders in time: EulerImplicit 1, BDF2 2
add_test(NAME CheckAllenCahn2DPeriodictest3_bdf2TimeConvergence
         COMMAND ${Python3_EXECUTABLE} time_convergence.py --check
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(CheckAllenCahn2DPeriodictest3_bdf2TimeConvergence PROPERTIES
                     DEPENDS AllenCahn2DPeriodictest3_bdf2 LABELS "2D AC")
                     
create_col_comparison("CompareAllenCahn2DPeriodictest3_bdf2TimeConvergence" "convergence_output_ref.csv"
                      "time_convergence_output.csv" -1 relativeEpsilon 1e-6 FALSE
                      "CheckAllenCahn2DPeriodictest3_bdf2TimeConvergence" "2D AC")

# Too expensive for the GitHub Actions Linux runners: disabled there only
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND DEFINED ENV{GITHUB_ACTIONS})
  set_tests_properties(AllenCahn2DPeriodictest3_bdf2
                       CheckAllenCahn2DPeriodictest3_bdf2TimeConvergence
                       CompareAllenCahn2DPeriodictest3_bdf2TimeConvergence
                       PROPERTIES DISABLED TRUE)
endif()
