# Create the library
option(SLOTH_USE_STATIC "Enable static library for Sloth" OFF)
if(SLOTH_USE_STATIC)
    add_library(Sloth STATIC)
else()
    add_library(Sloth SHARED)
endif(SLOTH_USE_STATIC)

file(GLOB_RECURSE SLOTH_SRC CONFIGURE_DEPENDS
    ${CMAKE_SOURCE_DIR}/kernel/src/*.cpp
)

target_sources(Sloth PRIVATE ${SLOTH_SRC})

set_target_properties(Sloth PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CMAKE_CXX_EXTENSIONS OFF
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
if(SLOTH_USE_EXTERNAL)
        target_link_directories(Sloth PUBLIC $ENV{EXT_LIBDIR} )
        target_link_libraries(Sloth PUBLIC $ENV{EXT_LIBNAME})
endif(SLOTH_USE_EXTERNAL)


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
install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION include
)

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/scripts/
    DESTINATION scripts
)
# SlothTargets.cmake creation
install(
    TARGETS Sloth
    EXPORT SlothTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include
)



install(
    EXPORT SlothTargets
    FILE SlothTargets.cmake
    NAMESPACE Sloth::
    DESTINATION cmake/Sloth
)

