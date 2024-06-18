create_test("Diffusion2Dtest1" "Diffusion2Dtest1" FALSE)

configure_file(${CMAKE_SOURCE_DIR}/kernel/Mesh-examples/star2D.msh ${CMAKE_CURRENT_BINARY_DIR}/star2D.msh COPYONLY)
