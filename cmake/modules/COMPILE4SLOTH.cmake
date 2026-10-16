
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

# Note: for the standard build types (Release, RelWithDebInfo, MinSizeRel),
# CMake already adds its own flags through CMAKE_CXX_FLAGS_<CONFIG>
# (Release: -O3 -DNDEBUG, RelWithDebInfo: -O2 -g -DNDEBUG,
#  MinSizeRel: -Os -DNDEBUG). They are not repeated below.

# Release
set(RELEASE_OPTIONS -DNO_BOUND_CHECK -march=native -flto)
if(NOT APPLE)
  # -fno-plt only makes sense for ELF binaries (Linux); no effect on macOS
  list(APPEND RELEASE_OPTIONS -fno-plt)
endif()

# Portable Release (custom build type: no CMAKE_CXX_FLAGS_PORTABLERELEASE,
# so optimisation flags are set here)
set(PORTABLE_REL_OPTIONS -O3 -DNDEBUG)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|i.86)$")
  # Generic x86-64 target
  list(APPEND PORTABLE_REL_OPTIONS -march=x86-64 -mtune=generic)
endif()
# On arm64/aarch64, no -march: the compiler default target is already generic
# (-march=x86-64 is invalid there)

# RelWithDebInfo
set(RELEASEWITHDEBINFO_OPTIONS -march=native -fno-omit-frame-pointer)

# MinSizeRel
set(MINSIZEREL_OPTIONS -flto)

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
        # --coverage at link time pulls the coverage runtime (libgcov with GCC,
        # profile runtime with Clang): no explicit -lgcov (absent with AppleClang)
        set_target_properties(${CURRENT_EXE} PROPERTIES LINK_FLAGS ${COVERAGE_LINK_OPTIONS})
    endif()
endfunction()
