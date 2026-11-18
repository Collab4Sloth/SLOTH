# ##############################################################################
# CUDA (MFEM built with CUDA, MFEM_USE_CUDA in config.mk)
# ##############################################################################
# The loops run on the device (mfem::forall with MFEM_HOST_DEVICE lambdas, in
# the headers of SLOTH: coefficients, integrators...) must be compiled by nvcc,
# as MFEM itself: the C++ sources of SLOTH and of the programs using it are then
# compiled as CUDA sources (see sloth_cuda_target, called by set_compile_options).
#
#   GPU architecture: -DCMAKE_CUDA_ARCHITECTURES=<XX> (e.g. 80 for A100), by
#   default the one of MFEM (-arch=sm_XX in MFEM_CXXFLAGS of config.mk).
#   Host compiler of nvcc: -DCMAKE_CUDA_HOST_COMPILER=<compiler>, by default the
#   C++ compiler of SLOTH.
#
# Without CUDA, nothing is done.
# ##############################################################################
if(MFEM_USE_CUDA)
  if(CMAKE_VERSION VERSION_LESS 3.25)
    message(FATAL_ERROR "MFEM is built with CUDA: CMake 3.25 or newer is required to build SLOTH "
                        "(found ${CMAKE_VERSION}).")
  endif()

  if(NOT CMAKE_CUDA_HOST_COMPILER)
    set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}")
  endif()
  if(NOT CMAKE_CUDA_ARCHITECTURES)
    _slothmfem_expand(_sloth_mfem_cxxflags "${_mk_MFEM_CXXFLAGS}")
    if(_sloth_mfem_cxxflags MATCHES "-arch[= ]+sm_([0-9]+)")
      set(CMAKE_CUDA_ARCHITECTURES "${CMAKE_MATCH_1}")
    else()
      message(FATAL_ERROR "MFEM is built with CUDA, but its GPU architecture was not found in "
                          "MFEM_CXXFLAGS (${MFEM_CONFIG_FILE}): use -DCMAKE_CUDA_ARCHITECTURES=<XX>.")
    endif()
  endif()

  enable_language(CUDA)
    
  # Run-time path of the CUDA libraries used by MFEM (cuSPARSE, cuBLAS...): the executables find
  # them without LD_LIBRARY_PATH (the libraries of the stubs directory are for linking only)
  find_package(CUDAToolkit REQUIRED)
  target_link_options(SlothMFEM INTERFACE "LINKER:-rpath,${CUDAToolkit_LIBRARY_DIR}")


  set(CMAKE_CUDA_STANDARD 20)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  set(CMAKE_CUDA_EXTENSIONS OFF)
  message(STATUS "MFEM: using CUDA, so SLOTH is compiled by ${CMAKE_CUDA_COMPILER} "
                 "(architecture ${CMAKE_CUDA_ARCHITECTURES}, host compiler ${CMAKE_CUDA_HOST_COMPILER})")

  # Options of nvcc required by the device loops of MFEM (as in MFEM's NVCC_FLAGS)
  set(SLOTH_CUDA_OPTIONS --expt-extended-lambda --expt-relaxed-constexpr)
endif()

# The C++ sources (.cpp) of a target are compiled as CUDA sources, with the options of nvcc.
# Called by set_compile_options for every target of SLOTH. No effect without CUDA.
function(sloth_cuda_target target)
  if(NOT MFEM_USE_CUDA)
    return()
  endif()
  get_target_property(_sources ${target} SOURCES)
  foreach(_source IN LISTS _sources)
    if(_source MATCHES "\\.cpp$")
      set_source_files_properties(${_source} TARGET_DIRECTORY ${target} PROPERTIES LANGUAGE CUDA)
    endif()
  endforeach()
  target_compile_options(${target} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:${SLOTH_CUDA_OPTIONS}>")
endfunction()

# Compile options of the host compiler: given as they are to the C++ compiler, through
# -Xcompiler to nvcc for the CUDA sources.
function(sloth_host_options target)
  if(NOT ARGN)
    return()
  endif()
  if(MFEM_USE_CUDA)
    list(JOIN ARGN "," _host)
    target_compile_options(${target} PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:${ARGN}>"
                                             "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${_host}>")
  else()
    target_compile_options(${target} PRIVATE ${ARGN})
  endif()
endfunction()
