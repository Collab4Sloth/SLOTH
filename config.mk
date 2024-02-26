# Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

# Variables corresponding to defines in config.hpp (YES, NO, or value)
MFEM_VERSION           = 40500
MFEM_VERSION_STRING    = 4.5
# see myEnv.sh
# MFEM_INSTALL_DIR       = 
MFEM_GIT_STRING        =
MFEM_USE_MPI           = YES
MFEM_USE_METIS         = YES
MFEM_USE_METIS_5       = YES
MFEM_DEBUG             = YES
MFEM_USE_EXCEPTIONS    = NO
MFEM_USE_ZLIB          = YES
MFEM_USE_LIBUNWIND     = NO
MFEM_USE_LAPACK        = NO
MFEM_THREAD_SAFE       = NO
MFEM_USE_LEGACY_OPENMP = NO
MFEM_USE_OPENMP        = YES
MFEM_USE_MEMALLOC      = YES
MFEM_TIMER_TYPE        = 2
MFEM_USE_SUNDIALS      = NO
MFEM_USE_MESQUITE      = NO
MFEM_USE_SUITESPARSE   = YES
MFEM_USE_SUPERLU       = YES
MFEM_USE_SUPERLU5      = $(MFEM_USE_SUPERLU5)
MFEM_USE_MUMPS         = NO
MFEM_USE_STRUMPACK     = YES
MFEM_USE_GINKGO        = NO
MFEM_USE_AMGX          = NO
MFEM_USE_GNUTLS        = NO
MFEM_USE_NETCDF        = NO
MFEM_USE_PETSC         = NO
MFEM_USE_SLEPC         = NO
MFEM_USE_MPFR          = NO
MFEM_USE_SIDRE         = NO
MFEM_USE_FMS           = NO
MFEM_USE_CONDUIT       = NO
MFEM_USE_PUMI          = NO
MFEM_USE_HIOP          = NO
MFEM_USE_GSLIB         = NO
MFEM_USE_CUDA          = NO
MFEM_USE_HIP           = NO
MFEM_USE_RAJA          = NO
MFEM_USE_OCCA          = NO
MFEM_USE_CEED          = NO
MFEM_USE_CALIPER       = NO
MFEM_USE_UMPIRE        = NO
MFEM_USE_SIMD          = NO
MFEM_USE_ADIOS2        = NO
MFEM_USE_MKL_CPARDISO  = NO
MFEM_USE_MOONOLITH     = NO
MFEM_USE_ADFORWARD     = NO
MFEM_USE_CODIPACK      = NO
MFEM_USE_BENCHMARK     = NO
MFEM_USE_PARELAG       = NO
MFEM_USE_ENZYME        = NO

# Compiler, compile options, and link options
MFEM_CXX       = $(OPENMPI)/bin/mpic++
MFEM_HOST_CXX  = $(MFEM_CXX)
MFEM_CPPFLAGS  =  -I$(SRC_DIR) 
MFEM_CXXFLAGS  = -g -Wall -std=c++17
MFEM_PUGUIXML  = #-I$(PUGIXML_INSTALL_DIR)/include
MFEM_TPLFLAGS  = \
				-I$(HYPRE)/include \
				-I$(SUPERLUDIST)/include \
				-I$(PARMETIS)/include \
				-I$(METIS)/include \
				-I$(SUNDIALS)/include \
				-I$(SUITE)/include \
				-I$(STRUM)/include -fopenmp \
				-I$(PARMETIS)/include \
				-I$(NETLIB)/include \
				-I$(BUTTERFLY)/include \
				-I$(ZFP)/include \
				-I$(PETSC)/include -fopenmp \
				-I$(ZLIB)/include
MFEM_INCFLAGS  = -I$(MFEM_INC_DIR) $(MFEM_TPLFLAGS) $(MFEM_PUGUIXML)
MFEM_PICFLAG   =
MFEM_FLAGS     = $(MFEM_CPPFLAGS) $(MFEM_CXXFLAGS) $(MFEM_INCFLAGS)
MFEM_EXT_LIBS  = -Wl,-rpath,$(HYPRE)/lib -Wl,-rpath,$(SPACK)/openblas-0.3.21-vyi6qw5caknyfrujfjgqatebpvzkm5kh/lib -L$(HYPRE)/lib -L$(SPACK)/openblas-0.3.21-vyi6qw5caknyfrujfjgqatebpvzkm5kh/lib -lHYPRE -lopenblas -Wl,-rpath,$(SUPERLUDIST)/lib -Wl,-rpath,$(PARMETIS)/lib -L$(SUPERLUDIST)/lib -L$(PARMETIS)/lib -lsuperlu_dist -lparmetis -Wl,-rpath,$(SPACK)/openblas-0.3.21-vyi6qw5caknyfrujfjgqatebpvzkm5kh/lib -L$(SPACK)/openblas-0.3.21-vyi6qw5caknyfrujfjgqatebpvzkm5kh/lib -lopenblas -Wl,-rpath,$(METIS)/lib -L$(METIS)/lib -lmetis -Wl,-rpath,$(SUNDIALS)/lib -L$(SUNDIALS)/lib -lsundials_arkode -lsundials_cvodes -lsundials_nvecserial -lsundials_kinsol -lsundials_nvecparallel -lsundials_nvecmpiplusx -Wl,-rpath,$(SUITE)/lib -L$(SUITE)/lib -lklu -lbtf -lumfpack -lcholmod -lcolamd -lamd -lcamd -lccolamd -lsuitesparseconfig -Wl,-rpath,$(STRUM)/lib -L$(STRUM)/lib -lstrumpack -Wl,-rpath,$(PARMETIS)/lib -L$(PARMETIS)/lib -lparmetis -Wl,-rpath,$(NETLIB)/lib -L$(NETLIB)/lib -lscalapack -Wl,-rpath,$(BUTTERFLY)/lib -L$(BUTTERFLY)/lib -ldbutterflypack -lzbutterflypack -Wl,-rpath,$(ZFP)/lib -L$(ZFP)/lib -lzfp -Wl,-rpath,$(PETSC)/lib -L$(PETSC)/lib -lpetsc  -lrt -Wl,-rpath,$(ZLIB)/lib -L$(ZLIB)/lib -lz
MFEM_LIBS      =  -L$(MFEM_LIB_DIR) -lmfem  $(MFEM_EXT_LIBS)  -lstdc++fs # -L$(PUGIXML_INSTALL_DIR)/lib -lpugixml
MFEM_LIB_FILE  = $(MFEM_LIB_DIR)/libmfem.a
MFEM_STATIC    = YES
MFEM_SHARED    = NO
MFEM_BUILD_TAG = Linux pleiades1077 x86_64
MFEM_PREFIX    = $(MFEM_INSTALL_DIR)
MFEM_INC_DIR   = $(MFEM_INSTALL_DIR)/include
MFEM_LIB_DIR   = $(MFEM_INSTALL_DIR)/lib

# Location of test.mk
MFEM_TEST_MK = $(MFEM_INSTALL_DIR)/share/mfem/test.mk

# Command used to launch MPI jobs
MFEM_MPIEXEC    = mpirun
MFEM_MPIEXEC_NP = -np
MFEM_MPI_NP     = 4

# The NVCC compiler cannot link with -x=cu
MFEM_LINK_FLAGS := $(filter-out -x=cu -xhip, $(MFEM_FLAGS))

# Optional extra configuration
