/**
 * @file Boundary.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Boundary class used to build a boundary with a name, an index, a type and, if needed, a
 * value
 * @version 0.1
 * @date 2025-09-05
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

#include "BCs/Boundary.hpp"

#include <limits>
#include <map>
#include <string>
#include <vector>

#include "Spatial/Spatial.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new Boundary::Boundary object (null value prescribed by default)
 *
 * @param boundary_name Name of the boundary.
 * @param boundary_index Index of the boundary.
 * @param boundary_type Type of the boundary conditions (Dirichlet, Neumann, Periodic,
 * Robin)
 *
 */
Boundary::Boundary(const std::string& boundary_name, int boundary_index,
                   const std::string& boundary_type)
    : boundary_name_(boundary_name),
      boundary_type_(BoundaryConditionType::from(boundary_type)),
      boundary_index_(boundary_index) {
  switch (boundary_type_) {
    case BoundaryConditionType::Dirichlet:
    case BoundaryConditionType::Neumann:
    case BoundaryConditionType::Robin:
      this->is_periodic_boundary_ = false;
      break;
    case BoundaryConditionType::Periodic:
      this->is_periodic_boundary_ = true;
      break;
    default:
      mfem::mfem_error(
          "Boundary::Boundary(): only Dirichlet, Neumann, Periodic and Robin BoundaryConditionType "
          "are available");
  }
}

/**
 * @brief Construct a new Boundary:: Boundary object
 *
 * @param boundary_name Name of the boundary.
 * @param boundary_index Index of the boundary.
 * @param boundary_type Type of the boundary conditions (Dirichlet, Neumann, Periodic, Robin)
 * @param boundary_value Value of the boundary condition.
 */
Boundary::Boundary(const std::string& boundary_name, int boundary_index,
                   const std::string& boundary_type, double boundary_value)
    : boundary_name_(boundary_name),
      boundary_type_(BoundaryConditionType::from(boundary_type)),
      boundary_index_(boundary_index),
      boundary_value_(boundary_value) {
  this->is_periodic_boundary_ = false;
  if (boundary_type_ != BoundaryConditionType::Dirichlet) {
    mfem::mfem_error(
        "Boundary::Boundary(): only Dirichlet BoundaryConditionType can be associated with a "
        "constant value.\n"
        "Robin and non-homogeneous Neumann BoundaryConditionType are managed with Coefficients "
        "and BoundaryIntegrators.\n"
        "Periodic and Homogeneous Neumann BoundaryConditionType do not need a value.\n");
  }
}

/**
 * @brief Return the index associated to the boundary
 *
 * @return int The index of the boundary.
 */
int Boundary::get_boundary_index() const { return this->boundary_index_; }

/**
 * @brief Return the boundary type associated to the boundary
 *
 * @return BoundaryConditionType::value The type of the boundary.
 */
BoundaryConditionType::value Boundary::get_boundary_type() const { return this->boundary_type_; }

/**
 * @brief Flag to identify periodic boundary
 *
 * @return true
 * @return false
 */
bool Boundary::is_periodic_boundary() const { return this->is_periodic_boundary_; }

/**
 * @brief Return the double value prescribed on boundary
 *
 * @return double The value prescribded.
 */
double Boundary::get_boundary_value() const { return this->boundary_value_; }

/**
 * @brief Destroy the Boundary:: Boundary object
 *
 */
Boundary::~Boundary() {}
