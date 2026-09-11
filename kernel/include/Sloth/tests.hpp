/**
 * @file tests.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Usefull aliases for tests depending on the dimension of space
 * @version 0.1
 * @date 2025-01-22
 *
 * @copyright CEA (C) 2025
 *
 * This file is part of SLOTH.
 *
 * SLOTH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLOTH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once
#include <vector>

#include "Factory/ListFactory.hpp"

/**
 * @brief Usefull aliases for tests depending on the dimension of space
 *
 * @tparam DIM
 */
template <int DIM>
struct Test {
  /// @brief Finite Element collection (mfem object)
  using FECollection = mfem::H1_FECollection;
  /// @brief Variables object (SLOTH object)
  using VARS = Variables<FECollection, DIM>;
  /// @brief Variable object (SLOTH object)
  using VAR = Variable<FECollection, DIM>;
  /// @brief Paraview Collection object (mfem object)
  using PSTCollection = mfem::ParaViewDataCollection;
  /// @brief Post-processing (SLOTH object)
  using PST = PostProcessing<FECollection, PSTCollection, DIM>;
  /// @brief Spatial Discretization object (SLOTH object)
  using SPA = SpatialDiscretization<mfem::H1_FECollection, DIM>;
  /// @brief Boundary condition object (SLOTH object)
  using BCS = BoundaryConditions<FECollection, DIM>;
};

template <>
struct Test<1> {
  /// @brief Finite Element collection (mfem object)
  using FECollection = mfem::H1_FECollection;
  /// @brief Variables object (SLOTH object)
  using VARS = Variables<FECollection, 1>;
  /// @brief Variable object (SLOTH object)
  using VAR = Variable<FECollection, 1>;
  /// @brief Paraview Collection object (mfem object)
  using PSTCollection = mfem::ParaViewDataCollection;
  /// @brief Post-processing (SLOTH object)
  using PST = PostProcessing<FECollection, PSTCollection, 1>;
  /// @brief Spatial Discretization object (SLOTH object)
  using SPA = SpatialDiscretization<mfem::H1_FECollection, 1>;
  /// @brief Boundary condition object (SLOTH object)
  using BCS = BoundaryConditions<FECollection, 1>;
};
template <>
struct Test<2> {
  /// @brief Finite Element collection (mfem object)
  using FECollection = mfem::H1_FECollection;
  /// @brief Variables object (SLOTH object)
  using VARS = Variables<FECollection, 2>;
  /// @brief Variable object (SLOTH object)
  using VAR = Variable<FECollection, 2>;
  /// @brief Paraview Collection object (mfem object)
  using PSTCollection = mfem::ParaViewDataCollection;
  /// @brief Post-processing (SLOTH object)
  using PST = PostProcessing<FECollection, PSTCollection, 2>;
  /// @brief Spatial Discretization object (SLOTH object)
  using SPA = SpatialDiscretization<mfem::H1_FECollection, 2>;
  /// @brief Boundary condition object (SLOTH object)
  using BCS = BoundaryConditions<FECollection, 2>;
};
template <>
struct Test<3> {
  /// @brief Finite Element collection (mfem object)
  using FECollection = mfem::H1_FECollection;
  /// @brief Variables object (SLOTH object)
  using VARS = Variables<FECollection, 3>;
  /// @brief Variable object (SLOTH object)
  using VAR = Variable<FECollection, 3>;
  /// @brief Paraview Collection object (mfem object)
  using PSTCollection = mfem::ParaViewDataCollection;
  /// @brief Post-processing (SLOTH object)
  using PST = PostProcessing<FECollection, PSTCollection, 3>;
  /// @brief Spatial Discretization object (SLOTH object)
  using SPA = SpatialDiscretization<mfem::H1_FECollection, 3>;
  /// @brief Boundary condition object (SLOTH object)
  using BCS = BoundaryConditions<FECollection, 3>;
};

/**
 * @brief Define the DIM, FE/mesh/post-processing aliases, the
 *        set(Periodic)SpatialDiscretization factories, and the
 *        scheme-independent problem aliases (PB_MPI, PB_CALPHAD,
 *        PB_PROPERTY - none of them depend on OPERATOR_TYPE, so unlike
 *        OPE/PB they are defined once here rather than duplicated per
 *        scheme) shared by every scheme (Transient/Steady) for one spatial
 *        dimension.
 */
#define SLOTH_DEFINE_COMMON(DIM_VALUE)                                                           \
  constexpr int DIM = DIM_VALUE;                                                                 \
  using FECollection = Test<DIM>::FECollection;                                                  \
  using VARS = Test<DIM>::VARS;                                                                  \
  using VAR = Test<DIM>::VAR;                                                                    \
  using PST = Test<DIM>::PST;                                                                    \
  using SPA = Test<DIM>::SPA;                                                                    \
  using SPAS = std::vector<SPA*>;                                                                \
  using BCS = Test<DIM>::BCS;                                                                    \
  using PB_MPI = MPI_Problem<VARS, PST>;                                                         \
  template <class CALPHAD>                                                                       \
  using PB_CALPHAD = Calphad_Problem<CALPHAD, VARS, PST>;                                        \
  template <class PROPERTY>                                                                      \
  using PB_PROPERTY = Property_problem<PROPERTY, VARS, PST>;                                     \
  template <class... Args>                                                                       \
  SPAS setSpatialDiscretization(std::size_t N, Args&&... args) {                                 \
    return ::setSpatialDiscretization<FECollection, DIM>(N, false, std::forward<Args>(args)...); \
  }                                                                                              \
  template <class... Args>                                                                       \
  SPAS setPeriodicSpatialDiscretization(std::size_t N, Args&&... args) {                         \
    return ::setSpatialDiscretization<FECollection, DIM>(N, true, std::forward<Args>(args)...);  \
  }

/**
 * @brief Re-import the DIM-common aliases from COMMON_NS into the current
 *        namespace, one at a time (using-declarations, not `using
 *        namespace`), so qualified lookups such as SlothTransient2D::DIM
 *        stay valid. Each imported name still refers to the same entity as
 *        in COMMON_NS, so importing from several namespaces at once is
 *        never ambiguous.
 */
#define SLOTH_IMPORT_COMMON(COMMON_NS)       \
  using COMMON_NS::DIM;                      \
  using COMMON_NS::FECollection;             \
  using COMMON_NS::VARS;                     \
  using COMMON_NS::VAR;                      \
  using COMMON_NS::PST;                      \
  using COMMON_NS::SPA;                      \
  using COMMON_NS::SPAS;                     \
  using COMMON_NS::BCS;                      \
  using COMMON_NS::PB_MPI;                   \
  using COMMON_NS::PB_CALPHAD;               \
  using COMMON_NS::PB_PROPERTY;              \
  using COMMON_NS::setSpatialDiscretization; \
  using COMMON_NS::setPeriodicSpatialDiscretization;

/**
 * @brief Define only the PREFIX-prefixed, scheme-specific PDE aliases
 *        (TransientOPE/TransientPB or SteadyOPE/SteadyPB) - the only two
 *        that actually depend on OPERATOR_TYPE. Does not define the
 *        unprefixed OPE/PB names, so it can be invoked more than once in
 *        the same namespace (one call per scheme) without any redefinition
 *        clash. FECollection/VARS/PST/DIM must already be visible where
 *        this macro is invoked.
 */
#define SLOTH_DEFINE_SCHEME_ONLY(OPERATOR_TYPE, PREFIX) \
  using PREFIX##OPE = OPERATOR_TYPE<FECollection, DIM>; \
  using PREFIX##PB = Problem<PREFIX##OPE, VARS, PST>;

/**
 * @brief Import the DIM-common aliases (including PB_MPI/PB_CALPHAD/
 *        PB_PROPERTY) from COMMON_NS, define the usual
 *        unprefixed PDE aliases (OPE, PB) for backward-compatible
 *        single-scheme usage, and additionally define their
 *        PREFIX-prefixed counterparts via SLOTH_DEFINE_SCHEME_ONLY.
 */
#define SLOTH_DEFINE_SCHEME(COMMON_NS, OPERATOR_TYPE, PREFIX) \
  SLOTH_IMPORT_COMMON(COMMON_NS)                              \
  using OPE = OPERATOR_TYPE<FECollection, DIM>;               \
  using PB = Problem<OPE, VARS, PST>;                         \
  SLOTH_DEFINE_SCHEME_ONLY(OPERATOR_TYPE, PREFIX)

/**
 * @brief 1D aliases: DIM-common (incl. PB_MPI/PB_CALPHAD/
 *        PB_PROPERTY, each defined once) + both PDE schemes,
 *        already prefixed (TransientPB, SteadyPB, ...) - a single
 *        `using namespace Sloth1D` brings in everything needed for mixed
 *        transient/steady usage.
 */
namespace Sloth1D {
SLOTH_DEFINE_COMMON(1)
SLOTH_DEFINE_SCHEME_ONLY(TransientOperator, Transient)
SLOTH_DEFINE_SCHEME_ONLY(SteadyOperator, Steady)
}  // namespace Sloth1D

/**
 * @brief 2D aliases: DIM-common (incl. PB_MPI/PB_CALPHAD/
 *        PB_PROPERTY, each defined once) + both PDE schemes,
 *        already prefixed (TransientPB, SteadyPB, ...) - a single
 *        `using namespace Sloth2D` brings in everything needed for mixed
 *        transient/steady usage.
 */
namespace Sloth2D {
SLOTH_DEFINE_COMMON(2)
SLOTH_DEFINE_SCHEME_ONLY(TransientOperator, Transient)
SLOTH_DEFINE_SCHEME_ONLY(SteadyOperator, Steady)
}  // namespace Sloth2D

/**
 * @brief 3D aliases: DIM-common (incl. PB_MPI/PB_CALPHAD/
 *        PB_PROPERTY, each defined once) + both PDE schemes,
 *        already prefixed (TransientPB, SteadyPB, ...) - a single
 *        `using namespace Sloth3D` brings in everything needed for mixed
 *        transient/steady usage.
 */
namespace Sloth3D {
SLOTH_DEFINE_COMMON(3)
SLOTH_DEFINE_SCHEME_ONLY(TransientOperator, Transient)
SLOTH_DEFINE_SCHEME_ONLY(SteadyOperator, Steady)
}  // namespace Sloth3D

/**
 * @brief 1D transient aliases (OPE, PB, ... unprefixed, plus TransientOPE,
 *        TransientPB, ... prefixed) for single-scheme usage.
 */
namespace SlothTransient1D {
SLOTH_DEFINE_SCHEME(Sloth1D, TransientOperator, Transient)
}
/**
 * @brief 2D transient aliases (OPE, PB, ... unprefixed, plus TransientOPE,
 *        TransientPB, ... prefixed) for single-scheme usage.
 */
namespace SlothTransient2D {
SLOTH_DEFINE_SCHEME(Sloth2D, TransientOperator, Transient)
}
/**
 * @brief 3D transient aliases (OPE, PB, ... unprefixed, plus TransientOPE,
 *        TransientPB, ... prefixed) for single-scheme usage.
 */
namespace SlothTransient3D {
SLOTH_DEFINE_SCHEME(Sloth3D, TransientOperator, Transient)
}

/**
 * @brief 1D steady aliases (OPE, PB, ... unprefixed, plus SteadyOPE,
 *        SteadyPB, ... prefixed) for single-scheme usage.
 */
namespace SlothSteady1D {
SLOTH_DEFINE_SCHEME(Sloth1D, SteadyOperator, Steady)
}
/**
 * @brief 2D steady aliases (OPE, PB, ... unprefixed, plus SteadyOPE,
 *        SteadyPB, ... prefixed) for single-scheme usage.
 */
namespace SlothSteady2D {
SLOTH_DEFINE_SCHEME(Sloth2D, SteadyOperator, Steady)
}

/**
 * @brief 3D steady aliases (OPE, PB, ... unprefixed, plus SteadyOPE,
 *        SteadyPB, ... prefixed) for single-scheme usage.
 */
namespace SlothSteady3D {
SLOTH_DEFINE_SCHEME(Sloth3D, SteadyOperator, Steady)
}

#undef SLOTH_DEFINE_COMMON
#undef SLOTH_IMPORT_COMMON
#undef SLOTH_DEFINE_SCHEME_ONLY
#undef SLOTH_DEFINE_SCHEME