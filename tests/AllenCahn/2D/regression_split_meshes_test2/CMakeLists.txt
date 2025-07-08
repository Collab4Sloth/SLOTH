#create_test("SplitMeshAllenCahn2Dtest2" "SplitMeshAllenCahn2Dtest2" FALSE "EXE AC" "--oversubscribe 4")
create_test("SplitMeshAllenCahn2Dtest2" "SplitMeshAllenCahn2Dtest2" FALSE "EXE AC" 2)

file(GLOB CAMEMBERT2D
  "${CMAKE_SOURCE_DIR}/kernel/Mesh-examples/camembert2D.00*"
)
file( COPY ${CAMEMBERT2D} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

create_col_comparison("CompareSplitMeshAllenCahn2Dtest2" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "SplitMeshAllenCahn2Dtest2" "AC")

