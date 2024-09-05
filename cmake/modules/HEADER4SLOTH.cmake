
# Header files 
set(HEADER_FILE "${CMAKE_SOURCE_DIR}/kernel/sloth.hpp")
set(TEMP_FILE "${CMAKE_SOURCE_DIR}/kernel/sloth_tmp.hpp")

# List of hpp files
file(GLOB_RECURSE HEADER_FILES "${CMAKE_SOURCE_DIR}/kernel/*.hpp")

# Remove sloth.hpp from the list
list(REMOVE_ITEM HEADER_FILES "${HEADER_FILE}")
if(NOT SLOTH_USE_OC)
    file(GLOB_RECURSE OC_HEADER_FILE ${CMAKE_SOURCE_DIR}/kernel/Calphad/GibbsEnergyMinimizer/OpenCalphad/*.hpp)
    foreach(HEADER ${OC_HEADER_FILE})
        list(REMOVE_ITEM HEADER_FILES "${HEADER}")
    endforeach() 
endif()




# Built the new sloth.hpp file (temporary file)/**
set(NEW_HEADER_CONTENT "/**\n * @file sloth.hpp\n * @author ci230846  (clement.introini@cea.fr)\n * @brief List of all header files contained in the kernel directory\n * @version 0.1\n * @date 2024-08-06\n *\n * Copyright CEA (c) 2024\n *\n*/\n")
foreach(HEADER ${HEADER_FILES})
    string(REPLACE "${CMAKE_SOURCE_DIR}/kernel/" "" RELATIVE_HEADER ${HEADER})
    set(NEW_HEADER_CONTENT "${NEW_HEADER_CONTENT}#include \"${RELATIVE_HEADER}\"\n")
endforeach()
# set(NEW_HEADER_CONTENT "${NEW_HEADER_CONTENT}\n#pragma once \n")

file(WRITE ${TEMP_FILE} "${NEW_HEADER_CONTENT}")

# Read the current sloth.hpp file
if(EXISTS ${HEADER_FILE})
    file(READ ${HEADER_FILE} CURRENT_HEADER_CONTENT)
else()
    set(CURRENT_HEADER_CONTENT "")
endif()

# Read the new sloth.hpp file
file(READ ${TEMP_FILE} NEW_HEADER_CONTENT_TEMP)

# CComparison of the two sloth.hpp files
if(NOT "${CURRENT_HEADER_CONTENT}" STREQUAL "${NEW_HEADER_CONTENT_TEMP}")
    # Replace the current sloth.hpp by the new one
    file(RENAME ${TEMP_FILE} ${HEADER_FILE})
    message(STATUS "Updated ${HEADER_FILE}")
else()
    # Delete the temporary sloth.hpp if equal to the current one.
    file(REMOVE ${TEMP_FILE})
    message(STATUS "${HEADER_FILE} is already up to date")
endif()
