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

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Spatial/Spatial.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Check the boundary type and set the periodic flag
 *
 */
void Boundary::check_type() {
  switch (this->boundary_type_) {
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
 * @brief Construct a new Boundary object identified by the name of a boundary attribute set of
 * the mesh (null value prescribed by default)
 *
 * @param boundary_name Name of the boundary attribute set.
 * @param boundary_type Type of the boundary conditions (Dirichlet, Neumann, Periodic, Robin)
 */
Boundary::Boundary(const std::string& boundary_name, const std::string& boundary_type)
    : boundary_name_(boundary_name), boundary_type_(BoundaryConditionType::from(boundary_type)) {
  this->check_type();
}

/**
 * @brief Construct a new Boundary object identified by the name of a boundary attribute set of
 * the mesh, with a constant Dirichlet value
 *
 * @param boundary_name Name of the boundary attribute set.
 * @param boundary_type Type of the boundary conditions (Dirichlet only)
 * @param boundary_value Value of the boundary condition.
 */
Boundary::Boundary(const std::string& boundary_name, const std::string& boundary_type,
                   double boundary_value)
    : boundary_name_(boundary_name),
      boundary_type_(BoundaryConditionType::from(boundary_type)),
      boundary_value_(boundary_value) {
  if (this->boundary_type_ != BoundaryConditionType::Dirichlet) {
    mfem::mfem_error(
        "Boundary::Boundary(): only Dirichlet BoundaryConditionType can be associated with a "
        "constant value.\n"
        "Robin and non-homogeneous Neumann BoundaryConditionType are managed with Coefficients "
        "and BoundaryIntegrators.\n"
        "Periodic and Homogeneous Neumann BoundaryConditionType do not need a value.\n");
  }
}

/**
 * @brief Construct a new Boundary::Boundary object (null value prescribed by default)
 *
 * @param boundary_name Name of the boundary.
 * @param boundary_index Index of the boundary (attribute - 1).
 * @param boundary_type Type of the boundary conditions (Dirichlet, Neumann, Periodic,
 * Robin)
 *
 */
Boundary::Boundary(const std::string& boundary_name, int boundary_index,
                   const std::string& boundary_type)
    : Boundary(boundary_name, boundary_type) {
  this->boundary_index_ = boundary_index;
  this->has_user_index_ = true;
  this->boundary_indices_ = {boundary_index};
}

/**
 * @brief Construct a new Boundary:: Boundary object
 *
 * @param boundary_name Name of the boundary.
 * @param boundary_index Index of the boundary (attribute - 1).
 * @param boundary_type Type of the boundary conditions (Dirichlet only)
 * @param boundary_value Value of the boundary condition.
 */
Boundary::Boundary(const std::string& boundary_name, int boundary_index,
                   const std::string& boundary_type, double boundary_value)
    : Boundary(boundary_name, boundary_type, boundary_value) {
  this->boundary_index_ = boundary_index;
  this->has_user_index_ = true;
  this->boundary_indices_ = {boundary_index};
}

/**
 * @brief Compute the indices (attribute - 1) covered by the boundary from the boundary attribute
 * sets of the mesh.
 *
 * - a set named after the boundary exists: its attributes are used; if an index is also given, the
 *   set must contain exactly the attribute index + 1;
 * - no such set: the index given by the user is used; without index, the boundary is undefined.
 *
 * @param bdr_sets Boundary attribute sets of the mesh.
 */
void Boundary::resolve_indices(const mfem::AttributeSets& bdr_sets) {
  if (bdr_sets.AttributeSetExists(this->boundary_name_)) {
    const mfem::Array<int>& attrs = bdr_sets.GetAttributeSet(this->boundary_name_);
    if (attrs.Size() == 0) {
      const std::string msg = "Boundary::resolve_indices(): boundary attribute set '" +
                              this->boundary_name_ + "' is empty.";
      mfem::mfem_error(msg.c_str());
    }
    std::set<int> ids;
    for (const int a : attrs) ids.insert(a - 1);

    if (this->has_user_index_ && (ids.size() != 1 || *ids.begin() != this->boundary_index_)) {
      std::string list;
      for (const int a : attrs) list += " " + std::to_string(a);
      const std::string msg = "Boundary::resolve_indices(): boundary '" + this->boundary_name_ +
                              "' is given with index " + std::to_string(this->boundary_index_) +
                              " (attribute " + std::to_string(this->boundary_index_ + 1) +
                              ") but the boundary attribute set of the mesh with this name " +
                              "contains the attribute(s)" + list +
                              ". Remove the index or make it consistent with the mesh.";
      mfem::mfem_error(msg.c_str());
    }
    this->boundary_indices_.assign(ids.begin(), ids.end());
    return;
  }

  if (!this->has_user_index_) {
    std::string names;
    for (const auto& n : bdr_sets.GetAttributeSetNames()) names += " '" + n + "'";
    if (names.empty()) names = " none";
    const std::string msg = "Boundary::resolve_indices(): boundary '" + this->boundary_name_ +
                            "' is given without index and is not a boundary attribute set of the "
                            "mesh. Available boundary attribute sets:" +
                            names;
    mfem::mfem_error(msg.c_str());
  }
  this->boundary_indices_ = {this->boundary_index_};
}

/**
 * @brief Return the name of the boundary
 *
 * @return const std::string& The name of the boundary.
 */
const std::string& Boundary::get_boundary_name() const { return this->boundary_name_; }

/**
 * @brief Flag to identify a boundary defined with an index
 *
 * @return true if an index is given by the user
 */
bool Boundary::has_user_index() const { return this->has_user_index_; }

/**
 * @brief Return the index associated to the boundary
 *
 * @return int The index of the boundary (-1 if the boundary is identified by its name only).
 */
int Boundary::get_boundary_index() const { return this->boundary_index_; }

/**
 * @brief Return the indices (attribute - 1) covered by the boundary
 *
 * @return const std::vector<int>& Sorted indices. Empty for a boundary identified by its name
 * only, before resolve_indices() is called.
 */
const std::vector<int>& Boundary::get_boundary_indices() const { return this->boundary_indices_; }

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
