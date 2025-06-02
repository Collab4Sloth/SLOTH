create_test("Diffusion3Dtest1" "Diffusion3Dtest1" FALSE "2D EXE DIFFU" 2)
configure_file(${CMAKE_SOURCE_DIR}/kernel/Mesh-examples/camembert3D.msh ${CMAKE_CURRENT_BINARY_DIR}/camembert3D.msh COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/tests/Diffusion/3D/test1/database/binary_database.h5 ${CMAKE_CURRENT_BINARY_DIR}/binary_database.h5 COPYONLY)

