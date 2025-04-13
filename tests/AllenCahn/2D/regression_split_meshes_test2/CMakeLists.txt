#create_test("SplitMeshAllenCahn2Dtest2" "SplitMeshAllenCahn2Dtest2" FALSE "EXE AC" "--oversubscribe 4")

set(CURRENT_EXE "SplitMeshAllenCahn2Dtest2")
add_executable(${CURRENT_EXE} main.cpp)
target_link_libraries(${CURRENT_EXE} Sloth)
add_test(NAME "SplitMeshAllenCahn2Dtest2" COMMAND mpirun --oversubscribe -np 4 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE})

file(GLOB CAMEMBERT2D
  "${CMAKE_SOURCE_DIR}/kernel/Mesh-examples/camembert2D.*"
)
file( COPY ${CAMEMBERT2D} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

#add_test(NAME CompareSplitMeshAllenCahn2Dtest2 COMMAND mpirun -n 4 SplitMeshAllenCahn2Dtest2)
#create_col_comparison("CompareSplitMeshAllenCahn2Dtest2" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "SplitMeshAllenCahn2Dtest2" "AC")

