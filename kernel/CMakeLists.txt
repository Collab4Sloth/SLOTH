# Create the library
add_library(Sloth SHARED)


file(GLOB_RECURSE SLOTH_SRC CONFIGURE_DEPENDS
    ${CMAKE_SOURCE_DIR}/kernel/src/*.cpp
)

target_sources(Sloth PRIVATE ${SLOTH_SRC})

set_target_properties(Sloth PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
)

# And the most important part set the include directories
target_include_directories(Sloth
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/kernel/include>
        $<INSTALL_INTERFACE:include>
)


# Include MFEM + ExprTk
target_include_directories(Sloth SYSTEM PUBLIC ${MFEM_DIR}/include $ENV{ExprTk_INCLUDE_DIR})

target_link_directories(Sloth PUBLIC ${MFEM_DIR}/lib  $ENV{HYPRE_DIR}/lib $ENV{MPI_DIR}/lib $ENV{PETSC_DIR}/lib $ENV{METIS_DIR}/lib )

target_link_libraries(Sloth PUBLIC ${MFEM_LIBRARY} ${MPI_CXX_LIBRARIES} z ${HYPRE_LIBRARY} ${SuiteSparse_LIBRARIES} ${METIS_LIBRARY} ${FILESYSTEM_LIB})


# PETSC ??
if(MFEM_USE_PETSC)
    target_link_libraries(Sloth PUBLIC petsc )
endif(MFEM_USE_PETSC)


# EXTERNAL??
option(SLOTH_USE_EXTERNAL "Enable the use of an external package that requires linking with external library" OFF)
option(SLOTH_USE_AUTO_EXTERNAL "Enable the use of an autonomous external package" OFF)
if(SLOTH_USE_EXTERNAL OR SLOTH_USE_AUTO_EXTERNAL)
    file(GLOB_RECURSE EXT_SRC_FILES $ENV{EXT_SRC}/*.hpp)
    add_library(SlothExternal PUBLIC ${EXT_SRC_FILES})
    if(SLOTH_USE_EXTERNAL)
        target_link_directories(Sloth PUBLIC $ENV{EXT_LIBDIR} )
        target_link_libraries(Sloth PUBLIC SlothExternal $ENV{EXT_LIBNAME})
    endif(SLOTH_USE_EXTERNAL)
endif()


# LibTorch API ??
if(SLOTH_USE_LIBTORCH)
    find_package(Torch REQUIRED)
    MESSAGE(STATUS "Sloth use libTorch API ")
    target_include_directories(Sloth PUBLIC ${TORCH_INCLUDE_DIRS} )
    target_link_libraries(Sloth PUBLIC ${TORCH_LIBRARIES})
    target_compile_definitions(Sloth PUBLIC SLOTH_USE_LIBTORCH=1)
endif(SLOTH_USE_LIBTORCH)


# ExprTk ??
option(SLOTH_USE_EXPRTK "Sloth use ExprTk" OFF)

if(SLOTH_USE_EXPRTK)
    target_compile_definitions(Sloth PUBLIC SLOTH_USE_EXPRTK=1)
endif()

## Debug option
option(MATOOLS_VERBOSITY_LEVEL_1 "add verbosity level 1" OFF)

if(MATOOLS_VERBOSITY_LEVEL_1)
  target_compile_definitions(Sloth PUBLIC MATIMERS_VEROBSITY_LEVEL_1)
endif()


# ########################
# # Installation 
# ########################
# SlothTargets.cmake creation

install(
    TARGETS Sloth
    EXPORT SlothTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION include
)
install(
    EXPORT SlothTargets
    FILE SlothTargets.cmake
    NAMESPACE Sloth::
    DESTINATION lib/cmake/Sloth
)


install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/scripts/
    DESTINATION scripts
)


install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/tests/
    DESTINATION tests
)