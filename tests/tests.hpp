/**
 * @file tests.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Usefull aliases for tests depending on the dimension of space
 * @version 0.1
 * @date 2025-01-22
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

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
