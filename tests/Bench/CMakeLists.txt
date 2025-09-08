set(CURRENT_EXE kissing_bubbles)
add_executable(${CURRENT_EXE} main.cpp)
target_link_libraries(${CURRENT_EXE} Sloth)
add_test(NAME BenchKissingBubblesTest1 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} -r 0 -nx 16 -ny 16 -nz 16 --duration 0.05 --test-case 0)
## Duration limit with n = 32 and it does not converge with n = 16
#add_test(NAME BenchKissingBubblesTest2 COMMAND mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_EXE} -r 0 -nx 32 -ny 32 -nz 32 --duration 0.05 --test-case 1)
set_compile_options(${CURRENT_EXE})
install(TARGETS ${CURRENT_EXE})
