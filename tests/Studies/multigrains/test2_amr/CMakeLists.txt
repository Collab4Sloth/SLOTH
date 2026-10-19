execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
configure_file(${CMAKE_SOURCE_DIR}/tests/Studies/multigrains/test2_amr/initial-grain.txt ${CMAKE_CURRENT_BINARY_DIR}/initial-grain.txt COPYONLY)
create_test("AMRGrains3Dtest1" "AMRGrains3Dtest1" FALSE "2D EXE AC" 2)

