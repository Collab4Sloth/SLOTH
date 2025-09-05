/*
 * Copyright © CEA 2022
 *
 * \brief Variables class used to build and manage variables
 *
 * \file Variables.hpp
 * \author ci230846
 * \date 28/03/2023
 */
#include <any>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Variables/Variable.hpp"
#pragma once

/**
 * @brief Class used to manage a list of Variable
 *
 */
template <class T, int DIM>
class Variables {
 private:
  T* fecollection_;
  std::vector<Variable<T, DIM>> vect_variables_;

 public:
  template <class... Args>
  explicit Variables(Args... args);

  Variables();

  void add(Variable<T, DIM> var);
  std::vector<Variable<T, DIM>> getVariables() const;
  size_t get_variables_number() const;
  Variable<T, DIM>& getIVariable(const int& i);
  Variable<T, DIM>& get_variable(const std::string& name);
  std::map<std::string, mfem::ParGridFunction> get_map_gridfunction() const;
  std::map<std::string, Variable<T, DIM>> get_map_variable() const;
  ~Variables();
};

/**
 * @brief Construct a new Variables:: Variables object
 *
 */
template <class T, int DIM>
Variables<T, DIM>::Variables() {}

/**
 * @brief Construct a new Variables:: Variables object
 *
 * @tparam Args
 * @param args
 */
template <class T, int DIM>
template <class... Args>
Variables<T, DIM>::Variables(Args... args) {
  this->vect_variables_ = std::vector<Variable<T, DIM>>{args...};
}

/**
 * @brief Add a new variable
 *
 * @param var variable to add
 */
template <class T, int DIM>
void Variables<T, DIM>::add(Variable<T, DIM> var) {
  this->vect_variables_.emplace_back(var);
}

/**
 * @brief get vector of variables
 *
 * @return std::vector<Variable>
 */
template <class T, int DIM>
std::vector<Variable<T, DIM>> Variables<T, DIM>::getVariables() const {
  return this->vect_variables_;
}

/**
 * @brief Return the number of variables
 *
 * @tparam T
 * @tparam DIM
 * @return size_t
 */
template <class T, int DIM>
size_t Variables<T, DIM>::get_variables_number() const {
  return this->vect_variables_.size();
}

/**
 * @brief get i-th variable
 *
 * @param i
 * @return Variable&
 */
template <class T, int DIM>
Variable<T, DIM>& Variables<T, DIM>::getIVariable(const int& i) {
  return this->vect_variables_[i];
}

/**
 * @brief return a map of GridFunction for each variable name
 *
 * @return std::map<std::string, mfem::ParGridFunction>
 */
template <class T, int DIM>
std::map<std::string, mfem::ParGridFunction> Variables<T, DIM>::get_map_gridfunction() const {
  std::map<std::string, mfem::ParGridFunction> map_var;
  for (auto vv : this->vect_variables_) {
    const std::string& name = vv.getVariableName();
    auto gf = vv.get_gf();
    map_var.try_emplace(name, gf);
  }
  return map_var;
}
/**
 * @brief return a map of variables for each variable name
 *
 * @return std::map<std::string, Variable>
 */
template <class T, int DIM>
std::map<std::string, Variable<T, DIM>> Variables<T, DIM>::get_map_variable() const {
  std::map<std::string, Variable<T, DIM>> map_var;
  for (const auto& vv : this->vect_variables_) {
    const std::string& name = vv.getVariableName();
    map_var.try_emplace(name, vv);
  }
  return map_var;
}

/**
 * @brief return the variable called vname
 *
 * @param vname
 * @return Variable&
 */
template <class T, int DIM>
Variable<T, DIM>& Variables<T, DIM>::get_variable(const std::string& vname) {
  int id = 0;
  const auto vect_size = static_cast<int>(this->vect_variables_.size());
  for (auto i = 0; i < vect_size; i++) {
    const std::string& name = this->vect_variables_[i].getVariableName();
    if (name == vname) {
      id = i;
      break;
    }
  }
  return this->vect_variables_[id];
}

/**
 * @brief Destroy the Variables:: Variables object
 *
 */
template <class T, int DIM>
Variables<T, DIM>::~Variables() {}
