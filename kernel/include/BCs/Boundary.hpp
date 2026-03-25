/**
 * @file Boundary.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Boundary class used to build a boundary with a name, an index, a type and, if needed, a
 * value
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
 * @brief Class for defining a boundary with a name, an index, a type and, if needed, a value
 *
 */
class Boundary {
 private:
  std::string boundary_name_;
  std::string boundary_type_;
  int boundary_index_;
  bool is_essential_boundary_{false};
  bool is_periodic_boundary_{false};
  double boundary_value_{0.};

 public:
  Boundary(const std::string &boundary_name, const int boundary_index,
           const std::string &boundary_type);
  Boundary(const std::string &boundary_name, const int boundary_index,
           const std::string &boundary_type, const double &boundary_value);
  int get_boundary_index() const;
  bool is_essential_boundary() const;
  bool is_periodic_boundary() const;
  double get_boundary_value() const;

  ~Boundary();
};
