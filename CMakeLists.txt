cmake_minimum_required(VERSION 3.12)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON )

project(sloth VERSION 1.0 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS NO)

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
include_directories(${CMAKE_SOURCE_DIR})
add_subdirectory(kernel)

include(HEADER4SLOTH)


# ########################
# # Doc Configuration
# ########################
if(NOT TARGET doc)
  add_subdirectory(doc)
endif(NOT TARGET doc)

# ########################
# # Test Configuration
# ########################
include(COMPILE4SLOTH)



# ########################
# LibTorch API ??
# ########################
option(SLOTH_USE_LIBTORCH "Enable the use of libtorch API" OFF)
if(SLOTH_USE_LIBTORCH)
    find_package(Torch REQUIRED)
endif(SLOTH_USE_LIBTORCH)



include(TEST4SLOTH)

include(CTest)
enable_testing()

add_subdirectory(tests)

# ########################
# # Code Coverage
# ########################
IF(CMAKE_BUILD_TYPE MATCHES Coverage)
  # Collect CodeCoverage analysis with gcovr
  add_custom_target(cc
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/Coverage
    COMMAND gcovr -r ${CMAKE_SOURCE_DIR} --sort-uncovered --exclude-unreachable-branches --xml=${CMAKE_BINARY_DIR}/Coverage/coverage.xml --html --html-details -o ${CMAKE_BINARY_DIR}/Coverage/coverage.html
    COMMENT "Generate Html report for code coverage analysis (see Coverage/coverage.html)"
  )
ENDIF()

# ########################
# # Static analysis
# ########################
# Check if the lint target is already defined
if(NOT TARGET lint)
  find_program(CPP_LINT cpplint)

  # Add a custom target for linting
  add_custom_target(lint
    COMMAND ${CPP_LINT}
    --filter=-runtime/references,-build/header_guard,-runtime/string
    --counting=total
    --extensions=cpp,cxx,hpp,h,hxx
    --linelength=100
    --recursive ${CMAKE_SOURCE_DIR}/kernel ${CMAKE_SOURCE_DIR}/tests
    COMMENT "Linting source files"
  )
endif()
