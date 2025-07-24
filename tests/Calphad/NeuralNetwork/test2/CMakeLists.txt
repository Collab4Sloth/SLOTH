create_test("CalphadInformedNeuralNetworktest2" "CalphadInformedNeuralNetworktest2" FALSE "EXE CALPHAD" 2)
configure_file(${CMAKE_SOURCE_DIR}/kernel/Mesh-examples/camembert2D.msh ${CMAKE_CURRENT_BINARY_DIR}/camembert2D.msh COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/model.pt ${CMAKE_CURRENT_BINARY_DIR}/model.pt COPYONLY)
