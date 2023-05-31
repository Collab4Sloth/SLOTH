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
MFEM_SOURCE_DIR        = /tmp/ci230846/spack-stage/spack-stage-mfem-4.5.0-xjxq2muqa5emqgtk7nde63gl4rxvxc6h/spack-src
MFEM_INSTALL_DIR       = /home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/mfem-4.5.0-xjxq2muqa5emqgtk7nde63gl4rxvxc6h
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
MFEM_CXX       = /home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/openmpi-4.1.4-d6lksze6gp4ujqtvvomuylw63efqopie/bin/mpic++
MFEM_HOST_CXX  = /home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/openmpi-4.1.4-d6lksze6gp4ujqtvvomuylw63efqopie/bin/mpic++
MFEM_CPPFLAGS  =  -I$(SRC_DIR)
MFEM_CXXFLAGS  = -g -Wall -std=c++17
MFEM_PUGUIXML  = #-I$(PUGIXML_INSTALL_DIR)/include
MFEM_TPLFLAGS  = -I/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/hypre-2.27.0-ce52bddurkaaupytcxnmbibvuyr6dm2f/include -I/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/superlu-dist-8.1.2-e6uwf67etsnifcqfrsajfqhtjhcbp4py/include -I/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/parmetis-4.0.3-56eufzeiw4c7oka4fo5atqrpb7pegfmg/include -I/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/metis-5.1.0-wsod6lfx6henngwcw66imtodafykpl2c/include -I/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/sundials-6.4.1-wbkipfmjhc4i4ejjr3ocbkvnxlhwlbdu/include -I/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/suite-sparse-5.13.0-ujhxe6inqyn6zgyvptsreckx2oc3jqoa/include -I/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/strumpack-7.0.1-erk4bgbiypyvr5kjr2hooewwnncqi5au/include -fopenmp -I/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/parmetis-4.0.3-56eufzeiw4c7oka4fo5atqrpb7pegfmg/include -I/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/netlib-scalapack-2.2.0-zcjmy24bd2a2gar4p42hhjv5paumpeb6/include -I/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/butterflypack-2.2.2-h726msqzes7vqptvcse2seqep4pbxkjt/include -I/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/zfp-0.5.5-b26deeqfywuv5uxvde2soyq3sx2scj5c/include -I/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/petsc-3.18.3-wa5ktxrtg5ipxslsrkvqxc3sumfbsmjj/include -fopenmp -I/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/zlib-1.2.13-pxz5yozsajxrqkol7hmfkgshayjsvfew/include
MFEM_INCFLAGS  = -I$(MFEM_INC_DIR) $(MFEM_TPLFLAGS) $(MFEM_PUGUIXML)
MFEM_PICFLAG   =
MFEM_FLAGS     = $(MFEM_CPPFLAGS) $(MFEM_CXXFLAGS) $(MFEM_INCFLAGS)
MFEM_EXT_LIBS  = -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/hypre-2.27.0-ce52bddurkaaupytcxnmbibvuyr6dm2f/lib -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/openblas-0.3.21-vyi6qw5caknyfrujfjgqatebpvzkm5kh/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/hypre-2.27.0-ce52bddurkaaupytcxnmbibvuyr6dm2f/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/openblas-0.3.21-vyi6qw5caknyfrujfjgqatebpvzkm5kh/lib -lHYPRE -lopenblas -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/superlu-dist-8.1.2-e6uwf67etsnifcqfrsajfqhtjhcbp4py/lib -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/parmetis-4.0.3-56eufzeiw4c7oka4fo5atqrpb7pegfmg/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/superlu-dist-8.1.2-e6uwf67etsnifcqfrsajfqhtjhcbp4py/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/parmetis-4.0.3-56eufzeiw4c7oka4fo5atqrpb7pegfmg/lib -lsuperlu_dist -lparmetis -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/openblas-0.3.21-vyi6qw5caknyfrujfjgqatebpvzkm5kh/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/openblas-0.3.21-vyi6qw5caknyfrujfjgqatebpvzkm5kh/lib -lopenblas -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/metis-5.1.0-wsod6lfx6henngwcw66imtodafykpl2c/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/metis-5.1.0-wsod6lfx6henngwcw66imtodafykpl2c/lib -lmetis -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/sundials-6.4.1-wbkipfmjhc4i4ejjr3ocbkvnxlhwlbdu/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/sundials-6.4.1-wbkipfmjhc4i4ejjr3ocbkvnxlhwlbdu/lib -lsundials_arkode -lsundials_cvodes -lsundials_nvecserial -lsundials_kinsol -lsundials_nvecparallel -lsundials_nvecmpiplusx -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/suite-sparse-5.13.0-ujhxe6inqyn6zgyvptsreckx2oc3jqoa/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/suite-sparse-5.13.0-ujhxe6inqyn6zgyvptsreckx2oc3jqoa/lib -lklu -lbtf -lumfpack -lcholmod -lcolamd -lamd -lcamd -lccolamd -lsuitesparseconfig -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/strumpack-7.0.1-erk4bgbiypyvr5kjr2hooewwnncqi5au/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/strumpack-7.0.1-erk4bgbiypyvr5kjr2hooewwnncqi5au/lib -lstrumpack -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/parmetis-4.0.3-56eufzeiw4c7oka4fo5atqrpb7pegfmg/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/parmetis-4.0.3-56eufzeiw4c7oka4fo5atqrpb7pegfmg/lib -lparmetis -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/netlib-scalapack-2.2.0-zcjmy24bd2a2gar4p42hhjv5paumpeb6/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/netlib-scalapack-2.2.0-zcjmy24bd2a2gar4p42hhjv5paumpeb6/lib -lscalapack -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/butterflypack-2.2.2-h726msqzes7vqptvcse2seqep4pbxkjt/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/butterflypack-2.2.2-h726msqzes7vqptvcse2seqep4pbxkjt/lib -ldbutterflypack -lzbutterflypack -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/zfp-0.5.5-b26deeqfywuv5uxvde2soyq3sx2scj5c/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/zfp-0.5.5-b26deeqfywuv5uxvde2soyq3sx2scj5c/lib -lzfp -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/petsc-3.18.3-wa5ktxrtg5ipxslsrkvqxc3sumfbsmjj/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/petsc-3.18.3-wa5ktxrtg5ipxslsrkvqxc3sumfbsmjj/lib -lpetsc  -lrt -Wl,-rpath,/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/zlib-1.2.13-pxz5yozsajxrqkol7hmfkgshayjsvfew/lib -L/home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/zlib-1.2.13-pxz5yozsajxrqkol7hmfkgshayjsvfew/lib -lz
MFEM_LIBS      =  -L$(MFEM_LIB_DIR) -lmfem  $(MFEM_EXT_LIBS) # -L$(PUGIXML_INSTALL_DIR)/lib -lpugixml
MFEM_LIB_FILE  = $(MFEM_LIB_DIR)/libmfem.a
MFEM_STATIC    = YES
MFEM_SHARED    = NO
MFEM_BUILD_TAG = Linux pleiades1077 x86_64
MFEM_PREFIX    = /home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/mfem-4.5.0-xjxq2muqa5emqgtk7nde63gl4rxvxc6h
MFEM_INC_DIR   = /home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/mfem-4.5.0-xjxq2muqa5emqgtk7nde63gl4rxvxc6h/include
MFEM_LIB_DIR   = /home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/mfem-4.5.0-xjxq2muqa5emqgtk7nde63gl4rxvxc6h/lib

# Location of test.mk
MFEM_TEST_MK = /home/ci230846/home-local/MyGitProjects/spack/opt/spack/linux-debian10-skylake/gcc-8.3.0/mfem-4.5.0-xjxq2muqa5emqgtk7nde63gl4rxvxc6h/share/mfem/test.mk

# Command used to launch MPI jobs
MFEM_MPIEXEC    = mpirun
MFEM_MPIEXEC_NP = -np
MFEM_MPI_NP     = 4

# The NVCC compiler cannot link with -x=cu
MFEM_LINK_FLAGS := $(filter-out -x=cu -xhip, $(MFEM_FLAGS))

# Optional extra configuration
