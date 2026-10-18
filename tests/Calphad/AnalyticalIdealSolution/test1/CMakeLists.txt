execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficient.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

create_test("CalphadAnalyticalIdealSolutiontest1" "CalphadAnalyticalIdealSolutiontest1" FALSE "EXE CALPHAD" 1)
configure_file(${CMAKE_SOURCE_DIR}/Mesh-examples/camembert2D.msh ${CMAKE_CURRENT_BINARY_DIR}/camembert2D.msh COPYONLY)
