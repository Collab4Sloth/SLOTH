# 1D Allen-Cahn: equilibrium of a planar interface (hyperbolic tangent profile).
#
# The double-well potential phi^2 (1 - phi)^2 is described in coefficients.json;
# the SLOTH generator (Python + SymPy) turns it into the C++ class DoubleWell
# (AllenCahnCoefficients.hpp, written in the build directory, with its first and
# second derivatives). The header is regenerated whenever coefficients.json changes:
# edit the expression, rebuild, rerun.
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/AllenCahnCoefficients.hpp
  COMMAND ${Python3_EXECUTABLE} ${Sloth_COEFFICIENT_SCRIPT} -r -f ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/coefficients.json
  COMMENT "Generating AllenCahnCoefficients.hpp from coefficients.json")

add_executable(allen_cahn main.cpp ${CMAKE_CURRENT_BINARY_DIR}/AllenCahnCoefficients.hpp)
target_include_directories(allen_cahn PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(allen_cahn PRIVATE Sloth::Sloth)

add_test(NAME allen_cahn_run COMMAND ${MPIEXEC_EXECUTABLE} -n 1 $<TARGET_FILE:allen_cahn>)
set_tests_properties(allen_cahn_run PROPERTIES FIXTURES_SETUP allen_cahn_results)

# Equilibrium check: tanh profile kept and interfacial energy equal to sigma
# (Python standard library only)
add_test(NAME allen_cahn_verify COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/verify.py)
set_tests_properties(allen_cahn_verify PROPERTIES FIXTURES_REQUIRED allen_cahn_results)
