create_test("CalphadAnalyticalTwoPhasetest1" "CalphadAnalyticalTwoPhasetest1" FALSE "EXE CALPHAD" 1)
configure_file(${CMAKE_SOURCE_DIR}/tests/Calphad/AnalyticalTwoPhase/test1/1d_mesh_non_uniform_coarse.msh ${CMAKE_CURRENT_BINARY_DIR}/1d_mesh_non_uniform_coarse.msh COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/tests/Calphad/AnalyticalTwoPhase/test1/1d_mesh_non_uniform_extra_fine.msh ${CMAKE_CURRENT_BINARY_DIR}/1d_mesh_non_uniform_extra_fine.msh COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/tests/Calphad/AnalyticalTwoPhase/test1/1d_mesh_non_uniform_mega_fine.msh ${CMAKE_CURRENT_BINARY_DIR}/1d_mesh_non_uniform_mega_fine.msh COPYONLY)
