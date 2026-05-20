if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

########################
# User Information
########################
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Debug build.")
elseif(CMAKE_BUILD_TYPE STREQUAL "Optim")
    message("Optim build.")
elseif(CMAKE_BUILD_TYPE STREQUAL "Coverage")
    message("Coverage build.")
elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
    message("MinSizeRel build.")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    message("RelWithDebInfo build.")
elseif(CMAKE_BUILD_TYPE STREQUAL "PortableRelease")
    message("PortableRelease build.")
else()
    message("Release build.")
endif()

########################
# Compile options
########################


set(WARNINGS -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wnon-virtual-dtor -Wold-style-cast -Wnull-dereference -Wdouble-promotion)
# Debug 
set(DEBUG_OPTIONS -g -fstack-protector -O0 -gdwarf-4 -pedantic ${WARNINGS})

# Release
set(RELEASE_OPTIONS -O3 -DNDEBUG -DNO_BOUND_CHECK -march=native -flto -fno-plt)

# Portable Release
set(PORTABLE_REL_OPTIONS -O3 -DNDEBUG -march=x86-64 -mtune=generic )

# RelWithDebInfo
set(RELEASEWITHDEBINFO_OPTIONS -O2 -g -march=native -fno-omit-frame-pointer)

# MinSizeRel
set(MINSIZEREL_OPTIONS -Os -DNDEBUG -flto )

# Optim
set(OPTIM_OPTIONS -g -O2)
# Coverage
set(COVERAGE_OPTIONS -g -O0 --coverage)
set(COVERAGE_LINK_OPTIONS --coverage)

# Function called at each build 
function(set_compile_options CURRENT_EXE)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(${CURRENT_EXE} PRIVATE ${DEBUG_OPTIONS})

    elseif(CMAKE_BUILD_TYPE STREQUAL "Optim")
        target_compile_options(${CURRENT_EXE} PRIVATE ${OPTIM_OPTIONS})

    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(${CURRENT_EXE} PRIVATE ${RELEASE_OPTIONS})

    elseif(CMAKE_BUILD_TYPE STREQUAL "PortableRelease")
        target_compile_options(${CURRENT_EXE} PRIVATE ${PORTABLE_REL_OPTIONS})

    elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        target_compile_options(${CURRENT_EXE} PRIVATE ${RELEASEWITHDEBINFO_OPTIONS})

    elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
        target_compile_options(${CURRENT_EXE} PRIVATE ${MINSIZEREL_OPTIONS})

    elseif(CMAKE_BUILD_TYPE STREQUAL "Coverage")
        target_compile_options(${CURRENT_EXE} PRIVATE ${COVERAGE_OPTIONS})
        link_libraries(gcov)
        set_target_properties(${CURRENT_EXE} PROPERTIES LINK_FLAGS ${COVERAGE_LINK_OPTIONS})
    endif()
endfunction()
