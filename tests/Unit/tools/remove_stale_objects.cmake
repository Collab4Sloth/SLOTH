# Removes the object directories of unit-test executables that no longer exist.
#
#   cmake -D OBJECT_ROOT=<build>/tests/Unit/CMakeFiles -D KEEP=Unit_A,Unit_B \
#         -P remove_stale_objects.cmake
#
# A module renamed or merged into another one leaves its directory
# CMakeFiles/Unit_<Name>.dir in the build tree, with *.gcno files compiled
# from old versions of the sources: gcovr cannot merge them with the current
# ones. Only the Unit_*.dir directories absent from KEEP are removed.
string(REPLACE "," ";" KEEP "${KEEP}")
file(GLOB object_dirs LIST_DIRECTORIES true "${OBJECT_ROOT}/Unit_*.dir")
foreach(dir IN LISTS object_dirs)
  get_filename_component(name "${dir}" NAME)
  string(REGEX REPLACE "\\.dir$" "" target "${name}")
  if(NOT target IN_LIST KEEP)
    message(STATUS "unit_coverage: removing stale ${name}")
    file(REMOVE_RECURSE "${dir}")
  endif()
endforeach()
