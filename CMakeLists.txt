cmake_minimum_required(VERSION 3.8)
project(sloth VERSION 1.0 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/INSTALL)

# Check GCC version
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)

  if(GCC_VERSION VERSION_LESS 9)
    find_package(Filesystem REQUIRED)
    set(ENV{FILESYSTEM_VAR} std::filesystem)
  else()
    set(ENV{FILESYSTEM_VAR})
  endif()
endif()

# ########################
# # Global Configuration
# ########################
# install
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/INSTALL)

# ########################
# # MFEM Configuration
# ########################
include(MFEM4SLOTH)

# ########################
# # Sloth header's only
# # library (interface)
# ########################
add_subdirectory(kernel)

# include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# ########################
# # Doc Configuration
# ########################
if(NOT TARGET doc)
  add_custom_target(doc)
  add_subdirectory(doc)
endif(NOT TARGET doc)

# ########################
# # Test Configuration
# ########################
include(TEST4SLOTH)

include(CTest)
enable_testing()

add_subdirectory(tests)
