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
