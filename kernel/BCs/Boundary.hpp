/*
 * Copyright © CEA 2023
 *
 * \brief Boundary class used to build a boundary with a name, an index, a type and, if needed, a
 * value
 *
 * \file Boundary.hpp
 * \author ci230846
 * \date 23/03/2023
 */
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "Spatial/Spatial.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"
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
 * @brief Construct a new Boundary:: Boundary object (null value prescribed by default)
 *
 * @param boundary_name
 * @param boundary_index
 * @param boundary_type
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
      throw std::runtime_error(
          "Boundary::Boundary(): only Dirichlet, Neumann, Periodic and Robin BoundaryConditionType "
          "are available");
      break;
  }
}

/**
 * @brief Construct a new Boundary:: Boundary object
 *
 * @param boundary_name
 * @param boundary_index
 * @param boundary_type
 * @param boundary_value
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
      throw std::runtime_error(
          "Boundary::Boundary(): only Dirichlet, Neumann, Periodic and Robin BoundaryConditionType "
          "are available");
      break;
  }
}

/**
 * @brief return the index associated to the boundary
 *
 * @return int
 */
int Boundary::get_boundary_index() const { return this->boundary_index_; }

/**
 * @brief flag to identify essential boundary
 *
 * @return true
 * @return false
 */
bool Boundary::is_essential_boundary() const { return this->is_essential_boundary_; }

/**
 * @brief flag to identify periodic boundary
 *
 * @return true
 * @return false
 */
bool Boundary::is_periodic_boundary() const { return this->is_periodic_boundary_; }

/**
 * @brief return the double value prescribed on boundary
 *
 * @return double
 */
double Boundary::get_boundary_value() const { return this->boundary_value_; }

/**
 * @brief Destroy the Boundary:: Boundary object
 *
 */
Boundary::~Boundary() {}
