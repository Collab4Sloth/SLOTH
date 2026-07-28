/**
 * @file BoundaryConditions.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief BoundaryConditions class used to build and manage boundary conditions
 * @version 0.1
 * @date 2025-09-05
 *
 * @anchor bcs
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

#include <limits>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include "BCs/Boundary.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Class used to manage boundary conditions
 *
 */
template <class T, int DIM>
class BoundaryConditions {
 private:
  T* fecollection_;
  mfem::ParFiniteElementSpace* fespace_;
  mfem::Array<int> Dirichlet_bdr_;
  mfem::Array<double> Dirichlet_value_;
  mfem::Array<int> Neumann_bdr_;
  mfem::Array<int> Robin_bdr_;
  mfem::Array<int> ess_tdof_list_;
  std::initializer_list<Boundary> boundaries_;

 public:
  template <class... Args>
  BoundaryConditions(SpatialDiscretization<T, DIM>* spatial, const Args&... boundaries);
  void SetDirichletBoundaryConditions(mfem::Vector& u, Coefficients& coefficients, std::vector<mfem::Vector> auxvars_unk);
  mfem::Array<int> GetEssentialDofs();
  ~BoundaryConditions();
  mfem::Array<int> get_marker_array(const std::string& boundary_type);
};

#include "BCs/BoundaryConditions.tpp"
