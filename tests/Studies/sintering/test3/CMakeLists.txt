# execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
#     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
create_test("Sintering3Dtest3" "Sintering3Dtest3" FALSE "3D EXE CHAC" 1)
