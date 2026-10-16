# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
# Contributions made by CEA, France

from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack.package import *


class Sloth(CMakePackage):
    """SLOTH: multiphase-field multicomponent framework based on MFEM (CEA).

    Once installed, SLOTH is used from another CMake project with
        find_package(Sloth REQUIRED)
        target_link_libraries(<target> PRIVATE Sloth::Sloth)
    (example installed in <prefix>/share/sloth/examples/minimal).
    """

    homepage = "https://github.com/Collab4Sloth/SLOTH"
    git = "https://github.com/Collab4Sloth/SLOTH.git"

    version("master", branch="master")
    version("2.1.0", tag="v2.1.0", preferred=True)
    version("2.0.1", tag="v2.0.1")
    version("2.0.0", tag="v2.0.0")

    variant("shared", default=False, description="Build Sloth as a shared library")
    variant("mumps", default=False, description="Enable the MUMPS direct solver (through MFEM)")
    variant("petsc", default=False, description="Use an MFEM built with PETSc")
    variant(
        "native",
        default=False,
        description="Optimise for the CPU of the build machine (-march=native and "
        "link-time optimisation). Keep it off for an installation shared between "
        "machines or used with another compiler.",
    )

    depends_on("cmake@3.13:", type="build")
    # SymPy is required at configure time (generation of coefficient classes)
    depends_on("python@3:", type="build")
    depends_on("py-sympy", type="build")

    depends_on("mpi")
    # SLOTH only needs a parallel MFEM with SuiteSparse; MFEM's other options
    # (SUNDIALS, SuperLU, ...) are accepted but not required.
    depends_on("mfem@4.9: +mpi +suite-sparse")
    depends_on("hypre@3.1: +shared")

    depends_on("mfem+shared", when="+shared")
    depends_on("mfem+mumps", when="+mumps")
    depends_on("mumps", when="+mumps")
    depends_on("mfem+petsc", when="+petsc")
    depends_on("hypre+int64", when="+petsc")
    depends_on("petsc+int64", when="+petsc")

    def cmake_args(self):
        return [
            # MFEM and its third-party libraries are read from MFEM's config.mk
            self.define("MFEM_DIR", self.spec["mfem"].prefix),
            self.define_from_variant("SLOTH_USE_SHARED", "shared"),
            self.define_from_variant("SLOTH_NATIVE_ARCH", "native"),
            self.define_from_variant("SLOTH_ENABLE_LTO", "native"),
        ]

    def setup_build_environment(self, env):
        # Kept for versions <= 2.1.0, whose build system reads these variables
        # (not needed by the current one, harmless otherwise)
        env.set("HYPRE_DIR", self.spec["hypre"].prefix)
        env.set("USER_INSTALL_PATH", self.prefix)
        if self.spec.satisfies("+mumps"):
            env.set("MUMPS_DIR", self.spec["mumps"].prefix)
