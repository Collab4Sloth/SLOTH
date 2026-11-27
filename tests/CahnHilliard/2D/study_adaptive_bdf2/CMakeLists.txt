# Study (not a non-regression test): spinodal decomposition (PFHub BM1) with BDF2 and adaptive
# time-stepping, up to t = 1e4. Not built by default: uncomment add_subdirectory(study_adaptive_bdf2)
# in tests/CahnHilliard/2D/CMakeLists.txt, then run by hand:
#   mpirun -np 4 ./CahnHilliard2Dstudy_adaptive_bdf2
execute_process(COMMAND ${Python3_EXECUTABLE} ${BUILD_SCRIPT_DIR}/GenerateCoefficient.py -r  -f  ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set(CURRENT_EXE CahnHilliard2Dstudy_adaptive_bdf2)
add_executable(${CURRENT_EXE} main.cpp)
target_link_libraries(${CURRENT_EXE} Sloth)
target_include_directories(${CURRENT_EXE}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/tests>
    PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
set_compile_options(${CURRENT_EXE})
add_dependencies(tests ${CURRENT_EXE})
