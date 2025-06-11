/*
 * Copyright © CEA 2022
 *
 * \brief Variable class used to build and manage variables of model
 *
 * \file Variable.hpp
 * \author ci230846
 * \date 19/01/2022
 */
#include <algorithm>
#include <any>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "BCs/BoundaryConditions.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

template <class T, int DIM>
class Variable {
 private:
  T* fecollection_;
  BoundaryConditions<T, DIM> bcs_;
  std::string variable_name_;
  mfem::ParFiniteElementSpace* fespace_;

  std::vector<std::tuple<std::string, std::vector<double>>> coordinates_;
  // std::shared_ptr<AnalyticalFunctions<DIM>> ics_;
  std::map<int, mfem::Vector> map_of_unk_;
  int depth_ = 2;
  mfem::Vector unk_;
  mfem::ParGridFunction uh_;

  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> analytical_solution_{nullptr};

  std::vector<std::string> additional_variable_info_;

  mfem::Array<int> el_attr_;

  void setVariableDepth(const int& depth);

  void add_variable_info(const std::string& var);
  std::function<double(const mfem::Vector&, double)> buildAnalyticalFunction(
      const AnalyticalFunctions<DIM>& analytical_function);

  void setInitialCondition(const AnalyticalFunctions<DIM>& initial_condition_name);
  void setInitialCondition(const mfem::FunctionCoefficient& initial_condition_function);
  void setInitialCondition(const double& initial_condition_value);

  void setAnalyticalSolution(const AnalyticalFunctions<DIM>& analytical_solution_name);
  void setAnalyticalSolution(const mfem::FunctionCoefficient& analytical_solution_function);
  void saveBeforeUpdate();
  void set_attributes(SpatialDiscretization<T, DIM>* spatial,
                      const std::set<std::string>& attribute_names);

 public:
  /////////////////////////////
  // Without attributes names   //
  /////////////////////////////
  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const AnalyticalFunctions<DIM>& initial_condition_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const AnalyticalFunctions<DIM>& initial_condition_name,
           const AnalyticalFunctions<DIM>& analytical_solution_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const AnalyticalFunctions<DIM>& initial_condition_name,
           const mfem::FunctionCoefficient& analytical_solution_function);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const mfem::FunctionCoefficient& initial_condition_function);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const mfem::FunctionCoefficient& initial_condition_function,
           const AnalyticalFunctions<DIM>& analytical_solution_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const mfem::FunctionCoefficient& initial_condition_function,
           const mfem::FunctionCoefficient& analytical_solution_function);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const double& initial_condition_value);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const double& initial_condition_value,
           const AnalyticalFunctions<DIM>& analytical_solution_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const double& initial_condition_value,
           const mfem::FunctionCoefficient& analytical_solution_function);

  /////////////////////////////
  // With attributes names   //
  /////////////////////////////
  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const AnalyticalFunctions<DIM>& initial_condition_name,
           const std::set<std::string>& attribute_names);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const AnalyticalFunctions<DIM>& initial_condition_name,
           const std::set<std::string>& attribute_names,
           const AnalyticalFunctions<DIM>& analytical_solution_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const AnalyticalFunctions<DIM>& initial_condition_name,
           const std::set<std::string>& attribute_names,
           const mfem::FunctionCoefficient& analytical_solution_function);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const mfem::FunctionCoefficient& initial_condition_function,
           const std::set<std::string>& attribute_names);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const mfem::FunctionCoefficient& initial_condition_function,
           const std::set<std::string>& attribute_names,
           const AnalyticalFunctions<DIM>& analytical_solution_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const mfem::FunctionCoefficient& initial_condition_function,
           const std::set<std::string>& attribute_names,
           const mfem::FunctionCoefficient& analytical_solution_function);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const double& initial_condition_value, const std::set<std::string>& attribute_names);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const double& initial_condition_value, const std::set<std::string>& attribute_names,
           const AnalyticalFunctions<DIM>& analytical_solution_name);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const int& depth,
           const double& initial_condition_value, const std::set<std::string>& attribute_names,
           const mfem::FunctionCoefficient& analytical_solution_function);

  template <class... Args>
  void set_additional_information(Args&&... add_var_info);

  mfem::Vector get_last() const;
  mfem::Vector get_second_to_last() const;
  std::string getVariableName() const;
  std::vector<std::string> get_additional_variable_info() const;
  // std::shared_ptr<AnalyticalFunctions<DIM>> getInitialCondition();
  void update(const mfem::Vector& unk);
  mfem::Vector get_unknown() const;
  std::map<int, mfem::Vector> get_map_unknown();
  mfem::ParGridFunction get_gf() const;
  mfem::ParGridFunction get_igf() const;
  // mfem::ParGridFunction get_analytical_solution();
  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> get_analytical_solution();
  BoundaryConditions<T, DIM>* get_boundary_conditions();
  mfem::ParFiniteElementSpace* get_fespace();
  std::vector<std::tuple<std::string, std::vector<double>>> get_coordinates();

  ~Variable();
};

////////////////////////////////
// Without attributes names   //
////////////////////////////////
/**
 * @brief Construct a new Variable:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth, const AnalyticalFunctions<DIM>& initial_condition_name)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();
  this->setVariableDepth(depth);

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();

  this->setInitialCondition(initial_condition_name);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_name
 * @param analytical_solution_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth, const AnalyticalFunctions<DIM>& initial_condition_name,
                           const AnalyticalFunctions<DIM>& analytical_solution_name)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  // std::apply([dim, initial_condition_name, this]() {
  this->setInitialCondition(initial_condition_name);
  // });
  this->setVariableDepth(depth);
  // std::apply([dim, analytical_solution_name, this]() {
  this->setAnalyticalSolution(analytical_solution_name);
  // });

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_name
 * @param analytical_solution_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth, const AnalyticalFunctions<DIM>& initial_condition_name,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();
  this->uh_.SetSpace(fespace_);

  const auto dim = spatial->get_dimension();
  std::apply(
      [dim, initial_condition_name, this]() { this->setInitialCondition(initial_condition_name); });
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth,
                           const mfem::FunctionCoefficient& initial_condition_function)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_function
 * @param analytical_solution_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth,
                           const mfem::FunctionCoefficient& initial_condition_function,
                           const AnalyticalFunctions<DIM>& analytical_solution_name)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);
  std::apply([dim, analytical_solution_name, this]() {
    this->setAnalyticalSolution(analytical_solution_name);
  });

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_function
 * @param analytical_solution_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth,
                           const mfem::FunctionCoefficient& initial_condition_function,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_value
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth, const double& initial_condition_value)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();
  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);

  this->additional_variable_info_.resize(0);
}
/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_value
 * @param analytical_solution_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth, const double& initial_condition_value,
                           const AnalyticalFunctions<DIM>& analytical_solution_name)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);
  std::apply([dim, analytical_solution_name, this]() {
    this->setAnalyticalSolution(analytical_solution_name);
  });

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_value
 * @param analytical_solution_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth, const double& initial_condition_value,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);

  this->additional_variable_info_.resize(0);
}

/////////////////////////////
// With attributes names   //
/////////////////////////////
/**
 * @brief Construct a new Variable:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth, const AnalyticalFunctions<DIM>& initial_condition_name,
                           const std::set<std::string>& attribute_names)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();
  this->setVariableDepth(depth);

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  this->set_attributes(spatial, attribute_names);
  // this->setInitialCondition(initial_condition_name);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_name
 * @param analytical_solution_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth, const AnalyticalFunctions<DIM>& initial_condition_name,
                           const std::set<std::string>& attribute_names,
                           const AnalyticalFunctions<DIM>& analytical_solution_name)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();

  // std::apply([dim, initial_condition_name, this]() {

  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_name);
  // });
  this->setVariableDepth(depth);
  // std::apply([dim, analytical_solution_name, this]() {
  this->setAnalyticalSolution(analytical_solution_name);
  // });

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_name
 * @param analytical_solution_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth, const AnalyticalFunctions<DIM>& initial_condition_name,
                           const std::set<std::string>& attribute_names,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();

  this->uh_.SetSpace(fespace_);

  const auto dim = spatial->get_dimension();

  this->set_attributes(spatial, attribute_names);
  std::apply(
      [dim, initial_condition_name, this]() { this->setInitialCondition(initial_condition_name); });
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth,
                           const mfem::FunctionCoefficient& initial_condition_function,
                           const std::set<std::string>& attribute_names)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();

  this->uh_.SetSpace(fespace_);
  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_function
 * @param analytical_solution_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth,
                           const mfem::FunctionCoefficient& initial_condition_function,
                           const std::set<std::string>& attribute_names,
                           const AnalyticalFunctions<DIM>& analytical_solution_name)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();

  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);
  std::apply([dim, analytical_solution_name, this]() {
    this->setAnalyticalSolution(analytical_solution_name);
  });

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_function
 * @param analytical_solution_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth,
                           const mfem::FunctionCoefficient& initial_condition_function,
                           const std::set<std::string>& attribute_names,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();

  this->uh_.SetSpace(fespace_);

  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_value
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth, const double& initial_condition_value,
                           const std::set<std::string>& attribute_names)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();
  this->uh_.SetSpace(fespace_);

  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);

  this->additional_variable_info_.resize(0);
}
/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_value
 * @param analytical_solution_name
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth, const double& initial_condition_value,
                           const std::set<std::string>& attribute_names,
                           const AnalyticalFunctions<DIM>& analytical_solution_name)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);
  std::apply([dim, analytical_solution_name, this]() {
    this->setAnalyticalSolution(analytical_solution_name);
  });

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param depth
 * @param initial_condition_value
 * @param analytical_solution_function
 */
template <class T, int DIM>

Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const int& depth, const double& initial_condition_value,
                           const std::set<std::string>& attribute_names,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->coordinates_ = spatial->get_coordinates();

  this->uh_.SetSpace(fespace_);

  this->set_attributes(spatial, attribute_names);
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);

  this->additional_variable_info_.resize(0);
}

/**
 * @brief Associate additional information to the variable
 *
 * @tparam T
 * @tparam DIM
 * @tparam Args
 * @param add_var_info
 */
template <class T, int DIM>
template <class... Args>
void Variable<T, DIM>::set_additional_information(Args&&... add_var_info) {
  if constexpr (sizeof...(add_var_info) == 0) {
    this->additional_variable_info_.resize(0);
  } else {
    (add_variable_info(std::forward<Args>(add_var_info)), ...);
  }
}

/**
 * @brief build the function associated to initial_condition_name
 *
 * @param initial_condition_name
 * @return std::function<double(const mfem::Vector&, double)>
 */
template <class T, int DIM>
std::function<double(const mfem::Vector&, double)> Variable<T, DIM>::buildAnalyticalFunction(
    const AnalyticalFunctions<DIM>& analytical_function_name) {
  // this->ics_ = std::make_shared<AnalyticalFunctions<DIM>>();
  std::shared_ptr<AnalyticalFunctions<DIM>> ics(
      new AnalyticalFunctions<DIM>(analytical_function_name));
  return ics->getFunction();
}

//////////////////////////////////
// Set initial conditions
//////////////////////////////////

/**
 * @brief Define an initial condition on the basis of an analytical function defined by its name
 *
 * @param initial_condition_name
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition(const AnalyticalFunctions<DIM>& initial_condition_name) {
  auto icf = this->buildAnalyticalFunction(initial_condition_name);
  mfem::FunctionCoefficient ic_fc(icf);
  mfem::VectorArrayCoefficient vc(1);
  vc.Set(0, &ic_fc, false);
  this->uh_ = 0.;
  if (this->el_attr_.Size() > 0) {
    for (int i = 0; i < this->el_attr_.Size(); i++) {
      this->uh_.ProjectCoefficient(vc, this->el_attr_[i]);
    }
  } else {
    this->uh_.ProjectCoefficient(ic_fc);
  }
  this->uh_.GetTrueDofs(this->unk_);
}

/**
 * @brief Define an initial condition on the basis of a FunctionCoefficient
 *
 * @param initial_condition_function
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition(
    const mfem::FunctionCoefficient& initial_condition_function) {
  auto ic_fc = initial_condition_function;
  this->uh_ = 0.;
  if (this->el_attr_.Size() > 0) {
    mfem::VectorArrayCoefficient vc(1);
    vc.Set(0, &ic_fc, false);
    for (int i = 0; i < this->el_attr_.Size(); i++) {
      this->uh_.ProjectCoefficient(vc, this->el_attr_[i]);
    }
  } else {
    this->uh_.ProjectCoefficient(ic_fc);
  }
  this->uh_.GetTrueDofs(this->unk_);
}

/**
 * @brief Define an initial condition on the basis of a double value
 *
 * @param initial_condition_value
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition(const double& initial_condition_value) {
  mfem::ConstantCoefficient ic_fc(initial_condition_value);
  this->uh_ = 0.;
  if (this->el_attr_.Size() > 0) {
    mfem::VectorArrayCoefficient vc(1);
    // std::cout << el_attr_.Size() << std::endl;
    vc.Set(0, &ic_fc, false);
    for (int i = 0; i < this->el_attr_.Size(); i++) {
      this->uh_.ProjectCoefficient(vc, this->el_attr_[i]);
    }
  } else {
    this->uh_.ProjectCoefficient(ic_fc);
  }
  this->uh_.GetTrueDofs(this->unk_);
}

//////////////////////////////////
//////////////////////////////////
/**
 * @brief Define an analytical solution
 *
 * @param analytical_solution_name
 */
template <class T, int DIM>
void Variable<T, DIM>::setAnalyticalSolution(
    const AnalyticalFunctions<DIM>& analytical_solution_name) {
  this->analytical_solution_ = std::make_shared<std::function<double(const mfem::Vector&, double)>>(
      this->buildAnalyticalFunction(analytical_solution_name));
}

/**
 * @brief Define an analytical solution
 *
 * @param analytical_solution_function
 */
template <class T, int DIM>
void Variable<T, DIM>::setAnalyticalSolution(
    const mfem::FunctionCoefficient& analytical_solution_function) {
  this->analytical_solution_ = std::make_shared<std::function<double(const mfem::Vector&, double)>>(
      analytical_solution_function);
}

/**
 * @brief Return the second-to-last term of the saved variables
 *
 * @tparam T
 * @tparam DIM
 * @return mfem::Vector
 */
template <class T, int DIM>
mfem::Vector Variable<T, DIM>::get_second_to_last() const {
  return std::prev(std::prev(this->map_of_unk_.end()))->second;
}

/**
 * @brief Return the last term of the saved variables
 *
 * @tparam T
 * @tparam DIM
 * @return mfem::Vector
 */
template <class T, int DIM>
mfem::Vector Variable<T, DIM>::get_last() const {
  return std::prev(this->map_of_unk_.end())->second;
}

/**
 * @brief update the GridFunction on the basis of its associated unknown vector
 *
 */
template <class T, int DIM>
void Variable<T, DIM>::update(const mfem::Vector& unk) {
  this->saveBeforeUpdate();
  this->unk_ = unk;
  this->uh_.SetFromTrueDofs(this->unk_);

  auto current_solution = std::prev(this->map_of_unk_.end());
  current_solution->second = this->unk_;
}

/**
 * @brief Save previous solutions before update
 *
 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM>
void Variable<T, DIM>::saveBeforeUpdate() {
  auto begin = this->map_of_unk_.begin();
  auto end = this->map_of_unk_.end();
  ++begin;
  for (auto it = begin; it != end; ++it) {
    auto itm = (it->first) - 1;
    this->map_of_unk_.at(itm) = it->second;
  }
}

/**
 * @brief return the unkown vector
 *
 * @return mfem::Vector
 *
 */
template <class T, int DIM>
mfem::Vector Variable<T, DIM>::get_unknown() const {
  return this->unk_;
}

/**
 * @brief Return the additionnal information associated to variable
 *
 * @tparam T
 * @tparam DIM
 * @return std::vector<std::string>
 */
template <class T, int DIM>
std::vector<std::string> Variable<T, DIM>::get_additional_variable_info() const {
  return this->additional_variable_info_;
}

/**
 * @brief Return a map of unknows
 *
 * @tparam T
 * @tparam DIM
 * @return std::map<int, mfem::Vector>
 */
template <class T, int DIM>
std::map<int, mfem::Vector> Variable<T, DIM>::get_map_unknown() {
  return this->map_of_unk_;
}

/**
 * @brief return the gridfunction associated to the unknown
 *
 * @return mfem::ParGridFunction
 */
template <class T, int DIM>
mfem::ParGridFunction Variable<T, DIM>::get_gf() const {
  return this->uh_;
}

/**
 * @brief return the function associated to the analytical solution
 *
 * @tparam T
 * @tparam DIM
 * @return std::function<double(const mfem::Vector&, double)>
 */
template <class T, int DIM>
std::shared_ptr<std::function<double(const mfem::Vector&, double)>>
Variable<T, DIM>::get_analytical_solution() {
  return this->analytical_solution_;
}

/**
 * @brief return the boundary condition object associated to the variable
 *
 * @tparam T
 * @tparam DIM
 * @return BoundaryConditions<T, DIM>
 */
template <class T, int DIM>
BoundaryConditions<T, DIM>* Variable<T, DIM>::get_boundary_conditions() {
  return &this->bcs_;
}

/**
 * @brief Get the name of the Variable
 *
 * @return std::string name of the variable
 */
template <class T, int DIM>
std::string Variable<T, DIM>::getVariableName() const {
  return this->variable_name_;
}

/**
 * @brief Set the variable depth and initialize at the given initial condition
 * @remark By default, 2 levels are considered. Not optimal in term of memory?
 * @remark By default, 3 levels in case of adaptive time-step?
 *
 * @tparam T
 * @tparam DIM
 * @param depth
 */
template <class T, int DIM>
void Variable<T, DIM>::setVariableDepth(const int& depth) {
  this->depth_ = std::max(2, depth);

  for (auto id = 0; id < depth; id++) {
    this->map_of_unk_.insert(std::pair<int, mfem::Vector>(id, this->unk_));
  }
}

/**
 * @brief
 *
 * @tparam T
 * @tparam DIM
 * @param var
 */
template <class T, int DIM>
void Variable<T, DIM>::add_variable_info(const std::string& var) {
  additional_variable_info_.emplace_back(var);
}

/**
 * @brief Return the pointer towards the FiniteElementSpace
 *
 */
template <class T, int DIM>
mfem::ParFiniteElementSpace* Variable<T, DIM>::get_fespace() {
  return this->fespace_;
}

/**
 * @brief Return the coordinates as {(X1,{x1,....,xn}),(X2,{y1,....,yn}),(X3,{z1,....,zn})}
 *
 * @tparam T
 * @tparam DIM
 * @return std::vector<std::tuple<std::string, std::vector<double>>>
 */
template <class T, int DIM>
std::vector<std::tuple<std::string, std::vector<double>>> Variable<T, DIM>::get_coordinates() {
  return this->coordinates_;
}

/**
 * @brief set the attributes for the variable
 *
 * @tparam T
 * @tparam DIM
 * @param spatial
 */
template <class T, int DIM>
void Variable<T, DIM>::set_attributes(SpatialDiscretization<T, DIM>* spatial,
                                      const std::set<std::string>& attribute_names) {
  std::shared_ptr<mfem::AttributeSets> elem_attr_sets = spatial->get_elem_attributes();
  std::set<std::string> set_mesh_attr_names = elem_attr_sets->GetAttributeSetNames();
  MFEM_VERIFY(std::includes(set_mesh_attr_names.begin(), set_mesh_attr_names.end(),
                            attribute_names.begin(), attribute_names.end()),
              "Error while setting attributes associated with elements for variable " +
                  this->variable_name_ + ". Please check your data");
  for (const auto& name : attribute_names) {
    this->el_attr_.Append(elem_attr_sets->GetAttributeSet(name));
  }
}

/**
 * @brief Destroy the Variable:: Variable object
 *
 */
template <class T, int DIM>
Variable<T, DIM>::~Variable() {}
