execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

create_test("SimpleCoupling2Dtest2_bdf2" "SimpleCoupling2Dtest2_bdf2" FALSE "2D EXE COUPLING" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/tools/time_convergence.py ${CMAKE_CURRENT_BINARY_DIR}/time_convergence.py COPYONLY)
# Time convergence of each problem of the partitioned coupling: implicit Euler order 1, BDF2
# order 2 (thanks to the extrapolation of the lagged auxiliary variable)
add_test(NAME CheckSimpleCoupling2Dtest2_bdf2TimeConvergence
         COMMAND bash -c "\"${Python3_EXECUTABLE}\" time_convergence.py -p AllenCahn -o time_convergence_eta --check && \"${Python3_EXECUTABLE}\" time_convergence.py -p CahnHilliard -o time_convergence_c --check"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(CheckSimpleCoupling2Dtest2_bdf2TimeConvergence PROPERTIES
                     DEPENDS SimpleCoupling2Dtest2_bdf2 LABELS "2D COUPLING")
# Non-regression on the error values 
create_col_comparison("CompareSimpleCoupling2Dtest2_bdf2TimeConvergenceEta" "time_convergence_eta.csv"
                      "time_convergence_eta.csv" -1 relativeEpsilon 1e-6 FALSE
                      "CheckSimpleCoupling2Dtest2_bdf2TimeConvergence" "2D COUPLING")
create_col_comparison("CompareSimpleCoupling2Dtest2_bdf2TimeConvergenceC" "time_convergence_c.csv"
                      "time_convergence_c.csv" -1 relativeEpsilon 1e-6 FALSE
                      "CheckSimpleCoupling2Dtest2_bdf2TimeConvergence" "2D COUPLING")
