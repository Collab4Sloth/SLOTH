create_test("CalphadInformedNeuralNetworktest4" "CalphadInformedNeuralNetworktest4" FALSE "EXE CALPHAD" 2)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/solid_model.pt ${CMAKE_CURRENT_BINARY_DIR}/solid_model.pt COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/liquid_model.pt ${CMAKE_CURRENT_BINARY_DIR}/liquid_model.pt COPYONLY)
