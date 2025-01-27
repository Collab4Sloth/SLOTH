if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

########################
# User Information
########################
if(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug build.")
elseif(CMAKE_BUILD_TYPE MATCHES Optim)
    message("Optim build.")
elseif(CMAKE_BUILD_TYPE MATCHES Coverage)
    message("Coverage build.")
else()
    message("Release build.")
endif()

########################
# Compile options
########################
set(WARNINGS -Wall -Wextra -Wno-deprecated -Wparentheses -Wreturn-type -Wmissing-declarations -fmessage-length=0 -Wunused -Wfatal-errors -Wpointer-arith -Wcast-align -Wwrite-strings -Wctor-dtor-privacy -Wnon-virtual-dtor -Woverloaded-virtual -Wfloat-equal -Wno-endif-labels -Wsign-compare -Wmissing-format-attribute -Wno-multichar -Wno-deprecated-declarations -Wpacked -Wredundant-decls -Wdisabled-optimization -Wunknown-pragmas -Wundef -Wreorder)
# Debug 
set(DEBUG_OPTIONS -g -O0 -gdwarf-4 -pedantic ${WARNINGS})
# Optim
set(OPTIM_OPTIONS -g -O2)
# Coverage
set(COVERAGE_OPTIONS -g -O0 --coverage)
set(COVERAGE_LINK_OPTIONS --coverage)

# Function called at each build 
function(set_compile_options target)
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        target_compile_options(${CURRENT_EXE} PRIVATE ${DEBUG_OPTIONS})

    elseif(CMAKE_BUILD_TYPE MATCHES Optim)
        target_compile_options(${CURRENT_EXE} PRIVATE ${OPTIM_OPTIONS})

    elseif(CMAKE_BUILD_TYPE MATCHES Coverage)
        target_compile_options(${CURRENT_EXE} PRIVATE ${COVERAGE_OPTIONS})
        link_libraries(gcov)
        set_target_properties(${CURRENT_EXE} PROPERTIES LINK_FLAGS ${COVERAGE_LINK_OPTIONS})
    endif()
endfunction()