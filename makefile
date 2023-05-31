# Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

# Use the MFEM build directory
# MFEM_DIR = /home/ci230846/home-local/MyGitProjects/MFEM-MGIS/mfem/
# MFEM_BUILD_DIR = /home/ci230846/home-local/MyGitProjects/MFEM-MGIS/mfem/build
# SRC = $($(MFEM_DIR,$(MFEM_DIR)/examples/,)
# CONFIG_MK = $(MFEM_BUILD_DIR)/config/config.mk
# Use the MFEM install directory
# MFEM_INSTALL_DIR = `spack location -i mfem`
# CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
CONFIG_MK = $(SRC_DIR)/config.mk

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

###################
# List of test
###################
## ALLEN-CAHN
ALLENCAHN_TESTS_1D= Tests/AllenCahn/1D/test1/main
ALLENCAHN_TESTS_2D= Tests/AllenCahn/2D/test1/main
ALLENCAHN_TESTS_2D_PERIODIC= Tests/AllenCahn/2D_periodic/test1/main\
							Tests/AllenCahn/2D_periodic/test2/main
ALLENCAHN_TESTS_3D= Tests/AllenCahn/3D/test1/main
ALLENCAHN_TESTS = $(ALLENCAHN_TESTS_1D) $(ALLENCAHN_TESTS_2D) $(ALLENCAHN_TESTS_2D_PERIODIC) $(ALLENCAHN_TESTS_3D) 
# ## CAHN-HILLIARD
# CAHNHILLIARD_TESTS_2D= Tests/CahnHilliard/2D/main
# CAHNHILLIARD_TESTS = $(CAHNHILLIARD_TESTS_2D) 

PF_TESTS = $(ALLENCAHN_TESTS) 
# $(CAHNHILLIARD)

###################
###################

SEQ_EXAMPLES = $(PF_TESTS)
PAR_EXAMPLES = $(PF_TESTS)

ifeq ($(MFEM_USE_MPI),NO)
   EXAMPLES = $(SEQ_EXAMPLES)
else
   EXAMPLES = $(PAR_EXAMPLES) $(SEQ_EXAMPLES)
endif
SUBDIRS =
ifeq ($(MFEM_USE_SUNDIALS),YES)
   SUBDIRS += sundials
endif
ifeq ($(MFEM_USE_PETSC),YES)
   SUBDIRS += petsc
endif
ifeq ($(MFEM_USE_PUMI),YES)
   SUBDIRS += pumi
endif
ifeq ($(MFEM_USE_HIOP),YES)
   SUBDIRS += hiop
endif
ifeq ($(MFEM_USE_GINKGO),YES)
   SUBDIRS += ginkgo
endif

SUBDIRS_ALL = $(addsuffix /all,$(SUBDIRS))
SUBDIRS_TEST = $(addsuffix /test,$(SUBDIRS))
SUBDIRS_CLEAN = $(addsuffix /clean,$(SUBDIRS))
SUBDIRS_TPRINT = $(addsuffix /test-print,$(SUBDIRS))

.SUFFIXES:
.SUFFIXES: .o *.cpp *.hpp *.h .mk
.PHONY: all clean clean-build clean-exec

# Remove built-in rule
%: %.cpp

# Replace the default implicit rule for *.cpp files
%: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

all:  $(EXAMPLES) $(SUBDIRS_ALL)

.PHONY: $(SUBDIRS_ALL) $(SUBDIRS_TEST) $(SUBDIRS_CLEAN) $(SUBDIRS_TPRINT)
$(SUBDIRS_ALL) $(SUBDIRS_TEST) $(SUBDIRS_CLEAN):
	$(MAKE) -C $(@D) $(@F)
$(SUBDIRS_TPRINT):
	@$(MAKE) -C $(@D) $(@F)

MFEM_TESTS = EXAMPLES
include $(MFEM_TEST_MK)
test: $(SUBDIRS_TEST)
test-print: $(SUBDIRS_TPRINT)


doc: 
	doxygen PhaseFieldDoxyfile

clean-doc:
	rm -fr docinfo latex html Doxyfile.bak doxygen.log 


# Testing: Parallel vs. serial runs
RUN_MPI = $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $(MFEM_MPI_NP)
%-test-par: %
	@$(call mfem-test,$<, $(RUN_MPI), Parallel example)
%-test-seq: %
	@$(call mfem-test,$<,, Serial example)

# Testing: "test" target and mfem-test* variables are defined in config/test.mk

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: clean-build clean-exec $(SUBDIRS_CLEAN)

clean-build:
	rm -f *.o *~ $(SEQ_EXAMPLES) $(PAR_EXAMPLES)
	rm -rf *.dSYM *.TVD.*breakpoints Output *.out *log

clean-exec:
	@rm -f main
