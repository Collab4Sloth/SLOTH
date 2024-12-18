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

#include "BCs/BoundaryConditions.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/AnalyticalFunctions.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

template <class T, int DIM>
class Variable {
 private:
  T* fecollection_;
  BoundaryConditions<T, DIM> bcs_;
  std::string variable_name_;
  VariableType::value variable_type_;
  mfem::ParFiniteElementSpace* fespace_;
  // std::shared_ptr<AnalyticalFunctions<DIM>> ics_;
  std::map<int, mfem::Vector> map_of_unk_;
  int depth_ = 2;
  mfem::Vector unk_;
  mfem::ParGridFunction uh_;

  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> analytical_solution_{nullptr};

  void setVariableDepth(const int& depth);

  std::function<double(const mfem::Vector&, double)> buildAnalyticalFunction(
      const int& dim, const AnalyticalFunctions<DIM>& analytical_function);

  void setInitialCondition(const int& dim, const AnalyticalFunctions<DIM>& initial_condition_name);
  void setInitialCondition(const mfem::FunctionCoefficient& initial_condition_function);
  void setInitialCondition(const double& initial_condition_value);

  void setAnalyticalSolution(const int& dim,
                             const AnalyticalFunctions<DIM>& analytical_solution_name);
  void setAnalyticalSolution(const mfem::FunctionCoefficient& analytical_solution_function);
  void saveBeforeUpdate();

 public:
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

  mfem::Vector get_last();
  std::string getVariableName() const;
  VariableType::value getVariableType();
  // std::shared_ptr<AnalyticalFunctions<DIM>> getInitialCondition();
  void update(const mfem::Vector& unk);
  mfem::Vector get_unknown();
  std::map<int, mfem::Vector> get_map_unknown();
  mfem::ParGridFunction get_gf() const;
  mfem::ParGridFunction get_igf() const;
  // mfem::ParGridFunction get_analytical_solution();
  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> get_analytical_solution();
  BoundaryConditions<T, DIM>* get_boundary_conditions();
  mfem::ParFiniteElementSpace* get_fespace();

  ~Variable();
};

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
  this->setVariableDepth(depth);

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();

  Variable<T, DIM>::setInitialCondition(dim, initial_condition_name);
  // mfem::ConstantCoefficient cc(0.);
  // this->uh_.ProjectCoefficient(cc);
  // this->uh_.GetTrueDofs(this->unk_);
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

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  // std::apply([dim, initial_condition_name, this]() {
  this->setInitialCondition(dim, initial_condition_name);
  // });
  this->setVariableDepth(depth);
  // std::apply([dim, analytical_solution_name, this]() {
  this->setAnalyticalSolution(dim, analytical_solution_name);
  // });
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

  this->uh_.SetSpace(fespace_);

  const auto dim = spatial->get_dimension();
  std::apply([dim, initial_condition_name, this]() {
    this->setInitialCondition(dim, initial_condition_name);
  });
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);
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

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);
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

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);
  std::apply([dim, analytical_solution_name, this]() {
    this->setAnalyticalSolution(dim, analytical_solution_name);
  });
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

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_function);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);
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
  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);
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

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);
  std::apply([dim, analytical_solution_name, this]() {
    this->setAnalyticalSolution(dim, analytical_solution_name);
  });
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

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_value);
  this->setVariableDepth(depth);
  this->setAnalyticalSolution(analytical_solution_function);
}

/**
 * @brief build the function associated to initial_condition_name
 *
 * @param initial_condition_name
 * @return std::function<double(const mfem::Vector&, double)>
 */
template <class T, int DIM>
std::function<double(const mfem::Vector&, double)> Variable<T, DIM>::buildAnalyticalFunction(
    const int& dim, const AnalyticalFunctions<DIM>& analytical_function_name) {
  // this->ics_ = std::make_shared<AnalyticalFunctions<DIM>>();
  std::shared_ptr<AnalyticalFunctions<DIM>> ics(
      new AnalyticalFunctions<DIM>(analytical_function_name));
  return ics->getFunction();
}

/**
 * @brief Define an initial condition on the basis of an analytical function defined by its name
 *
 * @param initial_condition_name
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition(const int& dim,
                                           const AnalyticalFunctions<DIM>& initial_condition_name) {
  auto icf = this->buildAnalyticalFunction(dim, initial_condition_name);
  mfem::FunctionCoefficient ic_fc(icf);
  this->uh_.ProjectCoefficient(ic_fc);
  this->uh_.GetTrueDofs(this->unk_);
  // this->uh_.SetTrueVector();
}

/**
 * @brief Define an initial condition on the basis of a FunctionCoefficient
 *
 * @param initial_condition_function
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition(
    const mfem::FunctionCoefficient& initial_condition_function) {
  auto icf = initial_condition_function;
  this->uh_.ProjectCoefficient(icf);
  this->uh_.GetTrueDofs(this->unk_);
  // this->uh_.SetTrueVector();
}

/**
 * @brief Define an initial condition on the basis of a double value
 *
 * @param initial_condition_value
 */
template <class T, int DIM>
void Variable<T, DIM>::setInitialCondition(const double& initial_condition_value) {
  mfem::ConstantCoefficient ic_fc(initial_condition_value);
  this->uh_.ProjectCoefficient(ic_fc);
  this->uh_.GetTrueDofs(this->unk_);
  // this->uh_.SetTrueVector();
}

/**
 * @brief Define an analytical solution
 *
 * @param analytical_solution_name
 */
template <class T, int DIM>
void Variable<T, DIM>::setAnalyticalSolution(
    const int& dim, const AnalyticalFunctions<DIM>& analytical_solution_name) {
  this->analytical_solution_ = std::make_shared<std::function<double(const mfem::Vector&, double)>>(
      this->buildAnalyticalFunction(dim, analytical_solution_name));
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
 * @brief Return the last term of the saved variables
 *
 * @tparam T
 * @tparam DIM
 * @return mfem::Vector
 */
template <class T, int DIM>
mfem::Vector Variable<T, DIM>::get_last() {
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
  // for (const auto& [k, v] : this->map_of_unk_) {
  //   std::cout << " value at it = " << k << std::endl;
  //   v.Print();
  // }
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
    this->map_of_unk_[itm] = it->second;
  }
}

/**
 * @brief return the unkown vector
 *
 * @return mfem::Vector
 *
 */
template <class T, int DIM>
mfem::Vector Variable<T, DIM>::get_unknown() {
  return this->unk_;
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
 * @brief return the gridfunction associated to the analytical solution
 *
 * @return mfem::ParGridFunction
 */
// template <class T, int DIM>
// mfem::ParGridFunction Variable<T, DIM>::get_analytical_solution() {
//   return this->uh_ex_;
// }

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
 * @brief Get the Type of the Variable
 *
 * @return VariableType::value type of the variable
 */
template <class T, int DIM>
VariableType::value Variable<T, DIM>::getVariableType() {
  return this->variable_type_;
}

// template <class T, int DIM>
// std::shared_ptr<const AnalyticalFunctions<DIM>&> Variable<T, DIM>::getInitialCondition() {
//   return std::shared_ptr<const AnalyticalFunctions<DIM>&>(this->ics_);
// }

/**
 * @brief Set the variable depth and initialize at the given initial condition
 * @remark By default, 2 levels are considered. Not optimal in term of memory?
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

// /**
//  * @brief Set the VariableType from string
//  *
//  * @param depth string defining the type of the variable
//  */
// template <class T, int DIM>
// void Variable<T, DIM>::setVariableType(const int& depth) {
//   switch (VariableType::from(type)) {
//     case VariableType::Conserved:
//       this->variable_type_ = VariableType::Conserved;
//       break;
//     case VariableType::Unconserved:
//       this->variable_type_ = VariableType::Unconserved;
//       break;
//     default:
//       mfem::mfem_error(
//           "Variable::getVariableType() : only conserved and unconserved type are allowed");
//       break;
//   }
// }

/**
 * @brief Return the pointer towards the FiniteElementSpace
 *
 */
template <class T, int DIM>
mfem::ParFiniteElementSpace* Variable<T, DIM>::get_fespace() {
  return this->fespace_;
}
/**
 * @brief Destroy the Variable:: Variable object
 *
 */
template <class T, int DIM>
Variable<T, DIM>::~Variable() {}
