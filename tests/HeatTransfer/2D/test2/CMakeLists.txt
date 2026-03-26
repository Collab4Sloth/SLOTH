execute_process(COMMAND python3 ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
create_test("HeatTransfer2Dtest2" "HeatTransfer2Dtest2" FALSE "2D EXE HEAT" 1)

add_test(NAME "VerifyHeatTransfer2Dtest2" COMMAND bash -c "python3 ${CMAKE_CURRENT_SOURCE_DIR}/verify.py" )
set_tests_properties(PROPERTIES WILL_FAIL FALSE LABELS "2D HEAT")
set_tests_properties(
    VerifyHeatTransfer2Dtest2
    PROPERTIES DEPENDS HeatTransfer2Dtest2
)

