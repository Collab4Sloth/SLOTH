/**
 * @file Boundary.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Boundary class used to build a boundary with a name, an index or an attribute set, a
 * type and, if needed, a value
 * @version 0.1
 * @date 2025-09-05
 *
 * Copyright CEA (C) 2025
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
#include <vector>

#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Class for defining a boundary with a name, a type and, if needed, a value.
 *
 * The boundary attributes it covers are given either by an index (attribute - 1), or by the name
 * of a boundary attribute set of the mesh (bdr_attribute_sets, MFEM >= 4.7). A set may gather
 * several attributes. When both an index and a set with the same name are available, they must
 * match.
 */
class Boundary {
 private:
  std::string boundary_name_;
  BoundaryConditionType::value boundary_type_;
  int boundary_index_{-1};
  bool has_user_index_{false};
  std::vector<int> boundary_indices_;
  bool is_periodic_boundary_{false};
  double boundary_value_{0.};

  void check_type();

 public:
  Boundary(const std::string& boundary_name, const std::string& boundary_type);
  Boundary(const std::string& boundary_name, const std::string& boundary_type,
           double boundary_value);
  Boundary(const std::string& boundary_name, int boundary_index, const std::string& boundary_type);
  Boundary(const std::string& boundary_name, int boundary_index, const std::string& boundary_type,
           double boundary_value);

  void resolve_indices(const mfem::AttributeSets& bdr_sets);

  const std::string& get_boundary_name() const;
  bool has_user_index() const;
  int get_boundary_index() const;
  const std::vector<int>& get_boundary_indices() const;
  BoundaryConditionType::value get_boundary_type() const;
  bool is_periodic_boundary() const;
  double get_boundary_value() const;

  ~Boundary();
};
