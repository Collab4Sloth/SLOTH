# ##############################################################################
# Tests of the Python scripts of SLOTH (scripts/), written with unittest
# (Python standard library: no dependency other than SymPy).
#
# Usage:
#   ctest -L Scripts
# or, without CMake:
#   cd tests/Scripts && PYTHONPATH=../../scripts python3 -m unittest -v
# ##############################################################################
add_test(NAME Scripts.VariationalForm
         COMMAND ${Python3_EXECUTABLE} -m unittest -v test_VariationalForm
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(Scripts.VariationalForm PROPERTIES
                     ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/scripts;PYTHONDONTWRITEBYTECODE=1"
                     LABELS "Scripts")
