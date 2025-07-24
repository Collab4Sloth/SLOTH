# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
# Contributions made by CEA, France

from spack import *
import os
import shutil

class Sloth(CMakePackage):
    """Sloth package"""

    homepage = "https://github.com/Collab4Sloth/SLOTH"
    url      = "https://github.com/Collab4Sloth/SLOTH.git"

    version('master', git='https://github.com/Collab4Sloth/SLOTH.git',  branch='master', preferred=True)
    version('1.0.0-alpha', git='https://github.com/Collab4Sloth/SLOTH.git',  tag='v1.0.0-alpha')
    version('1.0.0-alpha.1', git='https://github.com/Collab4Sloth/SLOTH.git',  tag='v1.0.0-alpha.1')

    variant('petsc'       , default=True  , description='Enable PETSc solvers, preconditioners, etc.')

    depends_on('hypre+int64', when='+petsc')
    depends_on('petsc+int64', when='+petsc')
    depends_on('mfem@4.7.0:+mpi+suite-sparse+sundials+superlu-dist+miniapps')
    depends_on('mfem@4.7.0:+petsc', when='+petsc')
    depends_on('cmake', type='build')

