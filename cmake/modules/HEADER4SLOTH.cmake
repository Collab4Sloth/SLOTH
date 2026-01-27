
# Header files 
set(HEADER_FILE "${CMAKE_SOURCE_DIR}/kernel/sloth.hpp")


# EXTERNAL ??
if(SLOTH_USE_EXTERNAL OR SLOTH_USE_AUTO_EXTERNAL)
    file(GLOB_RECURSE EXT_SRC_FILES $ENV{EXT_SRC}/*.hpp)
    foreach(HEADER ${EXT_SRC_FILES})
        set(NEW_HEADER_CONTENT "${NEW_HEADER_CONTENT}#include \"${HEADER}\"\n")
    endforeach()
endif()


file(APPEND ${HEADER_FILE} "${NEW_HEADER_CONTENT}")
