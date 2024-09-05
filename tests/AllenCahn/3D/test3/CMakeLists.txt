create_test("AllenCahn3Dtest3" "AllenCahn3Dtest3" FALSE "3D EXE AC" 8)

configure_file(${CMAKE_SOURCE_DIR}/kernel/Mesh-examples/camembert3D.msh ${CMAKE_CURRENT_BINARY_DIR}/camembert3D.msh COPYONLY)

