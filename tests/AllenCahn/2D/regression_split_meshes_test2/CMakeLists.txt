

file(GLOB CAMEMBERT2D
  "${CMAKE_SOURCE_DIR}/kernel/Mesh-examples/camembert2D.00*"
)
file( COPY ${CAMEMBERT2D} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)


set(CURRENT_EXE SplitMeshAllenCahn2Dtest2)
add_executable(${CURRENT_EXE} main.cpp)
target_link_libraries(${CURRENT_EXE} Sloth)

add_test(NAME SplitMeshAllenCahn2Dtest2 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${CURRENT_EXE}> ${MPIEXEC_POSTFLAGS})
set_tests_properties(SplitMeshAllenCahn2Dtest2 PROPERTIES WILL_FAIL FALSE LABELS "EXE AC" RUN_SERIAL TRUE)

set_compile_options(${CURRENT_EXE})
install(TARGETS ${CURRENT_EXE})

create_col_comparison("CompareSplitMeshAllenCahn2Dtest2" "time_specialized.csv" "Saves/Problem1/time_specialized.csv" -1 relativeEpsilon 1e-12 FALSE "SplitMeshAllenCahn2Dtest2" "AC")
