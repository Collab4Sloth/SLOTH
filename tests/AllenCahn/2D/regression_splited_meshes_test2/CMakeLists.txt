create_test("SplitedMeshAllenCahn2Dtest2" "SplitedMeshAllenCahn2Dtest2" FALSE "EXE AC" 3)


file(GLOB CAMEMBERT2D
  "${CMAKE_SOURCE_DIR}/kernel/Mesh-examples/camembert2D.*"
)
file( COPY ${CAMEMBERT2D} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

create_col_comparison("CompareSplitedMeshAllenCahn2Dtest2" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "SplitedMeshAllenCahn2Dtest2" "AC")

