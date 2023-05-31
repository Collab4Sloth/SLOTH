/*
 * Copyright © CEA 2022
 *
 * \brief Variable class used to build and manage variables of model
 *
 * \file Variable.hpp
 * \author ci230846
 * \date 19/01/2022
 */
#include <any>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include "BCs/BoundaryConditions.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/AnalyticalFunctions.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"
#pragma once

template <class T, int DIM>
class Variable {
 private:
  T* fecollection_;
  BoundaryConditions<T, DIM> bcs_;
  std::string variable_name_;
  VariableType::value variable_type_;
  mfem::FiniteElementSpace* fespace_;
  // std::shared_ptr<AnalyticalFunctions<DIM>> ics_;
  mfem::Vector unk_;
  mfem::GridFunction uh_;

  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> analytical_solution_{nullptr};

  void setVariableType(const std::string& type);
  template <class... Args>
  std::function<double(const mfem::Vector&, double)> buildAnalyticalFunction(
      const int& dim, const std::string& analytical_function, Args... args);

  template <class... Args>
  void setInitialCondition(const int& dim, const std::string& initial_condition_name, Args... args);
  void setInitialCondition(const mfem::FunctionCoefficient& initial_condition_function);
  void setInitialCondition(const double& initial_condition_value);

  template <class... Args>
  void setAnalyticalSolution(const int& dim, const std::string& analytical_solution_name,
                             Args... args);
  void setAnalyticalSolution(const mfem::FunctionCoefficient& analytical_solution_function);

 public:
  template <class... Args>
  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const std::string& type,
           const std::string& initial_condition_name, std::tuple<Args...> args1);

  template <class... Args1, class... Args2>
  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const std::string& type,
           const std::string& initial_condition_name, std::tuple<Args1...> args1,
           const std::string& analytical_solution_name, std::tuple<Args2...> args2);

  template <class... Args>
  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const std::string& type,
           const std::string& initial_condition_name, std::tuple<Args...> args1,
           const mfem::FunctionCoefficient& analytical_solution_function);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const std::string& type,
           const mfem::FunctionCoefficient& initial_condition_function);

  template <class... Args>
  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const std::string& type,
           const mfem::FunctionCoefficient& initial_condition_function,
           const std::string& analytical_solution_name, std::tuple<Args...> args1);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const std::string& type,
           const mfem::FunctionCoefficient& initial_condition_function,
           const mfem::FunctionCoefficient& analytical_solution_function);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const std::string& type,
           const double& initial_condition_value);

  template <class... Args>
  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const std::string& type,
           const double& initial_condition_value, const std::string& analytical_solution_name,
           std::tuple<Args...> args1);

  Variable(SpatialDiscretization<T, DIM>* spatial, const BoundaryConditions<T, DIM>& bcs,
           const std::string& variable_name, const std::string& type,
           const double& initial_condition_value,
           const mfem::FunctionCoefficient& analytical_solution_function);

  std::string getVariableName() const;
  VariableType::value getVariableType();
  // std::shared_ptr<AnalyticalFunctions<DIM>> getInitialCondition();
  void update(const mfem::Vector& unk);
  mfem::Vector get_unknown();
  mfem::GridFunction get_gf() const;
  mfem::GridFunction get_igf() const;
  // mfem::GridFunction get_analytical_solution();
  std::shared_ptr<std::function<double(const mfem::Vector&, double)>> get_analytical_solution();
  BoundaryConditions<T, DIM>* get_boundary_conditions();

  ~Variable();
};

/**
 * @brief Construct a new Variable:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param type
 * @param initial_condition_name
 */
template <class T, int DIM>
template <class... Args>
Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const std::string& type, const std::string& initial_condition_name,
                           std::tuple<Args...> args1)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->setVariableType(type);

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();

  std::apply(
      [dim, initial_condition_name, this](Args... args) {
        Variable<T, DIM>::setInitialCondition(dim, initial_condition_name, args...);
      },
      args1);
  // mfem::ConstantCoefficient cc(0.);
  // this->uh_.ProjectCoefficient(cc);
  // this->uh_.GetTrueDofs(this->unk_);
}

/**
 * @brief Construct a new Variable:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param type
 * @param initial_condition_name
 * @param analytical_solution_name
 */
template <class T, int DIM>
template <class... Args1, class... Args2>
Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const std::string& type, const std::string& initial_condition_name,
                           std::tuple<Args1...> args1, const std::string& analytical_solution_name,
                           std::tuple<Args2...> args2)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->setVariableType(type);

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  std::apply(
      [dim, initial_condition_name, this](Args1... args) {
        this->setInitialCondition(dim, initial_condition_name, args...);
      },
      args1);
  std::apply(
      [dim, analytical_solution_name, this](Args2... args) {
        this->setAnalyticalSolution(dim, analytical_solution_name, args...);
      },
      args2);
}

/**
 * @brief Construct a new Variable:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param type
 * @param initial_condition_name
 * @param analytical_solution_function
 */
template <class T, int DIM>
template <class... Args>
Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const std::string& type, const std::string& initial_condition_name,
                           std::tuple<Args...> args1,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->setVariableType(type);

  this->uh_.SetSpace(fespace_);

  const auto dim = spatial->get_dimension();
  std::apply([dim, initial_condition_name, this](
                 Args... args) { this->setInitialCondition(dim, initial_condition_name, args...); },
             args1);
  this->setAnalyticalSolution(analytical_solution_function);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param type
 * @param initial_condition_function
 */
template <class T, int DIM>
Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const std::string& type,
                           const mfem::FunctionCoefficient& initial_condition_function)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->setVariableType(type);

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_function);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param type
 * @param initial_condition_function
 * @param analytical_solution_name
 */
template <class T, int DIM>
template <class... Args>
Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const std::string& type,
                           const mfem::FunctionCoefficient& initial_condition_function,
                           const std::string& analytical_solution_name, std::tuple<Args...> args1)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->setVariableType(type);

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  this->setInitialCondition(initial_condition_function);
  std::apply(
      [dim, analytical_solution_name, this](Args... args) {
        this->setAnalyticalSolution(dim, analytical_solution_name, args...);
      },
      args1);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param type
 * @param initial_condition_function
 * @param analytical_solution_function
 */
template <class T, int DIM>
Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const std::string& type,
                           const mfem::FunctionCoefficient& initial_condition_function,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->setVariableType(type);

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_function);
  this->setAnalyticalSolution(analytical_solution_function);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param type
 * @param initial_condition_value
 */
template <class T, int DIM>
Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const std::string& type, const double& initial_condition_value)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->setVariableType(type);
  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_value);
}
/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param type
 * @param initial_condition_value
 * @param analytical_solution_name
 */
template <class T, int DIM>
template <class... Args>
Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const std::string& type, const double& initial_condition_value,
                           const std::string& analytical_solution_name, std::tuple<Args...> args1)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->setVariableType(type);

  this->uh_.SetSpace(fespace_);
  const auto dim = spatial->get_dimension();
  this->setInitialCondition(initial_condition_value);
  std::apply(
      [dim, analytical_solution_name, this](Args... args) {
        this->setAnalyticalSolution(dim, analytical_solution_name, args...);
      },
      args1);
}

/**
 * @brief Construct a new Variable<T>:: Variable object
 *
 * @param fespace
 * @param variable_name
 * @param type
 * @param initial_condition_value
 * @param analytical_solution_function
 */
template <class T, int DIM>
Variable<T, DIM>::Variable(SpatialDiscretization<T, DIM>* spatial,
                           const BoundaryConditions<T, DIM>& bcs, const std::string& variable_name,
                           const std::string& type, const double& initial_condition_value,
                           const mfem::FunctionCoefficient& analytical_solution_function)
    : bcs_(bcs), variable_name_(variable_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->setVariableType(type);

  this->uh_.SetSpace(fespace_);
  this->setInitialCondition(initial_condition_value);
  this->setAnalyticalSolution(analytical_solution_function);
}

/**
 * @brief build the function associated to initial_condition_name
 *
 * @param initial_condition_name
 * @return std::function<double(const mfem::Vector&, double)>
 */
template <class T, int DIM>
template <class... Args>
std::function<double(const mfem::Vector&, double)> Variable<T, DIM>::buildAnalyticalFunction(
    const int& dim, const std::string& analytical_function_name, Args... args) {
  // this->ics_ = std::make_shared<AnalyticalFunctions<DIM>>();
  std::shared_ptr<AnalyticalFunctions<DIM>> ics(new AnalyticalFunctions<DIM>);
  return ics->getAnalyticalFunctions(analytical_function_name, args...);
}

/**
 * @brief Define an initial condition on the basis of an analytical function defined by its name
 *
 * @param initial_condition_name
 */
template <class T, int DIM>
template <class... Args>
void Variable<T, DIM>::setInitialCondition(const int& dim,
                                           const std::string& initial_condition_name,
                                           Args... args) {
  auto icf = this->buildAnalyticalFunction(dim, initial_condition_name, args...);
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
template <class... Args>
void Variable<T, DIM>::setAnalyticalSolution(const int& dim,
                                             const std::string& analytical_solution_name,
                                             Args... args) {
  this->analytical_solution_ = std::make_shared<std::function<double(const mfem::Vector&, double)>>(
      this->buildAnalyticalFunction(dim, analytical_solution_name, args...));
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
 * @brief update the GridFunction on the basis of its associated unknown vector
 *
 */
template <class T, int DIM>
void Variable<T, DIM>::update(const mfem::Vector& unk) {
  this->unk_ = unk;
  this->uh_.SetFromTrueDofs(this->unk_);
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
 * @brief return the gridfunction associated to the unknown
 *
 * @return mfem::GridFunction
 */
template <class T, int DIM>
mfem::GridFunction Variable<T, DIM>::get_gf() const {
  return this->uh_;
}

/**
 * @brief return the gridfunction associated to the analytical solution
 *
 * @return mfem::GridFunction
 */
// template <class T, int DIM>
// mfem::GridFunction Variable<T, DIM>::get_analytical_solution() {
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
// std::shared_ptr<AnalyticalFunctions<DIM>> Variable<T, DIM>::getInitialCondition() {
//   return std::shared_ptr<AnalyticalFunctions<DIM>>(this->ics_);
// }

/**
 * @brief Set the VariableType from string
 *
 * @param type string defining the type of the variable
 */
template <class T, int DIM>
void Variable<T, DIM>::setVariableType(const std::string& type) {
  switch (VariableType::from(type)) {
    case VariableType::Conserved:
      this->variable_type_ = VariableType::Conserved;
      break;
    case VariableType::Unconserved:
      this->variable_type_ = VariableType::Unconserved;
      break;
    default:
      throw std::runtime_error(
          "Variable::getVariableType() : only conserved and unconserved type are allowed");
      break;
  }
}

/**
 * @brief Destroy the Variable:: Variable object
 *
 */
template <class T, int DIM>
Variable<T, DIM>::~Variable() {}
