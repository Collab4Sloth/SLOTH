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
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

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
  Boundary(const std::string &boundary_name, const int &boundary_index,
           const std::string &boundary_type);
  Boundary(const std::string &boundary_name, const int &boundary_index,
           const std::string &boundary_type, const double &boundary_value);
  int get_boundary_index() const;
  bool is_essential_boundary() const;
  bool is_periodic_boundary() const;
  double get_boundary_value() const;

  ~Boundary();
};

/**
 * @brief Construct a new Boundary::Boundary object (null value prescribed by default)
 *
 * @param boundary_name Name of the boundary.
 * @param boundary_index Index of the boundary.
 * @param boundary_type Type of the boundary conditions (Dirichlet, Neumann, Periodic,
 * Robin)
 *
 * @remark Robin BCs is not implemented yet
 */
Boundary::Boundary(const std::string &boundary_name, const int &boundary_index,
                   const std::string &boundary_type)
    : boundary_name_(boundary_name), boundary_index_(boundary_index) {
  switch (BoundaryConditionType::from(boundary_type)) {
    case BoundaryConditionType::Dirichlet:
      this->is_essential_boundary_ = true;
      this->is_periodic_boundary_ = false;
      break;
    case BoundaryConditionType::Neumann:
    case BoundaryConditionType::Robin:
      this->is_essential_boundary_ = false;
      this->is_periodic_boundary_ = false;
      break;
    case BoundaryConditionType::Periodic:
      this->is_periodic_boundary_ = true;
      this->is_essential_boundary_ = false;
      break;
    default:
      mfem::mfem_error(
          "Boundary::Boundary(): only Dirichlet, Neumann, Periodic and Robin BoundaryConditionType "
          "are available");
      break;
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
Boundary::Boundary(const std::string &boundary_name, const int &boundary_index,
                   const std::string &boundary_type, const double &boundary_value)
    : boundary_name_(boundary_name),
      boundary_index_(boundary_index),
      boundary_value_(boundary_value) {
  switch (BoundaryConditionType::from(boundary_type)) {
    case BoundaryConditionType::Dirichlet:
      this->is_essential_boundary_ = true;
      this->is_periodic_boundary_ = false;
      break;
    case BoundaryConditionType::Neumann:
    case BoundaryConditionType::Robin:
      this->is_essential_boundary_ = false;
      this->is_periodic_boundary_ = false;
      break;
    case BoundaryConditionType::Periodic:
      this->is_periodic_boundary_ = true;
      this->is_essential_boundary_ = false;
      break;
    default:
      mfem::mfem_error(
          "Boundary::Boundary(): only Dirichlet, Neumann, Periodic and Robin BoundaryConditionType "
          "are available");
      break;
  }
}

/**
 * @brief Return the index associated to the boundary
 *
 * @return int The index of the boundary.
 */
int Boundary::get_boundary_index() const { return this->boundary_index_; }

/**
 * @brief Flag to identify essential boundary
 *
 * @return true
 * @return false
 */
bool Boundary::is_essential_boundary() const { return this->is_essential_boundary_; }

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
