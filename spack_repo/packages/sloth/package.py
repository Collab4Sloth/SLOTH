# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
# Contributions made by CEA, France

from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack.package import *

import os
import shutil

class Sloth(CMakePackage):
    """Sloth package"""

    homepage = "https://github.com/Collab4Sloth/SLOTH"
    url      = "https://github.com/Collab4Sloth/SLOTH.git"

    version("master", git='https://github.com/Collab4Sloth/SLOTH.git', branch='master', preferred=True)
    version("1.0.0",  tag='v1.0.0')

    variant('petsc'       , default=False  , description='Enable PETSc solvers, preconditioners, etc.')
    variant('shared'       , default=False  , description='Enable Building MFEM with dynamic libraries.')

    depends_on('hypre+int64', when='+petsc')
    depends_on('petsc+int64', when='+petsc')
    depends_on('mfem@4.8.0:+mpi+suite-sparse+sundials+superlu-dist+miniapps')
    depends_on('mfem@4.8.0:+petsc', when='+petsc')
    depends_on('mfem@4.8.0:+shared', when='+shared')
    depends_on('cmake', type='build')

    def setup_build_environment(self, env):
        hypre_prefix = self.spec['hypre'].prefix
        env.set('HYPRE_DIR', hypre_prefix)
        env.set('USER_INSTALL_PATH', self.prefix)
