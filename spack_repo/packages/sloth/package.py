# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

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

    maintainers("CINTROINI")

    license("LGPL-3.0-or-later", checked_by="CINTROINI")

    version("master", branch="master")
    # Before the release: built from the Git tag. After tagging, replace with the
    # line given by "spack checksum sloth 2.2.0" (and add url = ... above)
    version("2.2.0", tag="v2.2.0")

    variant("shared", default=False, description="Build Sloth as a shared library")
    variant("mumps", default=False, description="Enable the MUMPS direct solver (through MFEM)")
    variant("voro", default=False, description="Enable Voro++ (Voronoi tessellations)")
    variant(
        "bigint",
        default=False,
        description="64-bit global indices (HYPRE +mixedint), for problems with more "
        "than about 2 billion unknowns",
    )
    variant(
        "native",
        default=False,
        description="Optimise for the CPU of the build machine (-march=native and "
        "link-time optimisation). Keep it off for an installation shared between "
        "machines or used with another compiler.",
    )

    depends_on("cxx", type="build")
    depends_on("cmake@3.13:", type="build")
    # SymPy: configure time (test coefficients) and run time (scripts/GenerateCoefficient.py)
    depends_on("python@3:", type=("build", "run"))
    depends_on("py-sympy", type=("build", "run"))

    depends_on("mpi")
    depends_on("mfem@4.9: +mpi")
    depends_on("hypre@3.1: +shared")

    depends_on("mfem+shared", when="+shared")
    depends_on("mfem+mumps", when="+mumps")
    depends_on("mumps", when="+mumps")
    depends_on("voropp", when="+voro")
    # mixedint rather than int64: 64-bit global indices, 32-bit local ones
    # (less memory, and compatible with HYPRE's external LAPACK)
    depends_on("hypre+mixedint", when="+bigint")

    def cmake_args(self):
        return [
            # MFEM and its third-party libraries are read from MFEM's config.mk
            self.define("MFEM_DIR", self.spec["mfem"].prefix),
            self.define_from_variant("SLOTH_USE_SHARED", "shared"),
            self.define_from_variant("SLOTH_USE_VORO", "voro"),
            self.define_from_variant("SLOTH_NATIVE_ARCH", "native"),
            self.define_from_variant("SLOTH_ENABLE_LTO", "native"),
        ]

    def setup_build_environment(self, env):
        # HYPRE_DIR: fallback used to read the HYPRE version (HYPRE4SLOTH.cmake).
        # USER_INSTALL_PATH: takes priority over CMAKE_INSTALL_PREFIX in
        # CMakeLists.txt; set to the Spack prefix so that a value exported in
        # the user's shell cannot redirect the installation.
        env.set("HYPRE_DIR", self.spec["hypre"].prefix)
        env.set("USER_INSTALL_PATH", self.prefix)
        # VORO: searched before /usr/local and Homebrew (kernel/CMakeLists.txt),
        # so that the Spack Voro++ is used even if another one is installed.
        if self.spec.satisfies("+voro"):
            env.set("VORO", self.spec["voropp"].prefix)

