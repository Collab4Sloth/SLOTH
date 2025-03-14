cmake_minimum_required(VERSION 3.22)

project(
  sloth-prototype-onika
  VERSION 0.0.1
  LANGUAGES CXX)

find_package(OpenMP REQUIRED)
find_package(onika REQUIRED)

add_subdirectory(src)

add_executable(SlothProto SlothProto.cpp)
target_compile_options(SlothProto PUBLIC ${ONIKA_COMPILE_OPTIONS})
target_compile_definitions(SlothProto PUBLIC ${ONIKA_COMPILE_DEFINITIONS})
target_include_directories(SlothProto PUBLIC ${ONIKA_INCLUDE_DIRS})
target_link_directories(SlothProto PUBLIC ${ONIKA_LIBRARY_DIRS})
target_link_libraries(SlothProto PUBLIC ${ONIKA_LIBRARIES})
