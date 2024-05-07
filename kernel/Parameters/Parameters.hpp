/*
 * Copyright © CEA 2022
 *
 * \brief Parameter class used to build and manage calculation parameter
 *
 * \file Parameter.hpp
 * \author ci230846
 * \date 19/01/2022
 */
#include <any>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "../Parameters/Parameter.hpp"
#pragma once

/**
 * @brief Class used to manage a list of Parameter
 *
 */
class Parameters {
 private:
  std::vector<Parameter> vect_params_;
  double get_val(const std::string& name) const;

 public:
  Parameters();
  template <class... Args>
  explicit Parameters(const Args&... args);

  void add(const Parameter& param);
  void ListParamByName();
  double get_parameter_value(const std::string& name) const;
  double get_parameter_value_or_default(const std::string& name, const double& default_value) const;
  bool get_option_value(const std::string& name) const;
  std::map<std::string, double> getMapParameters() const;
  ~Parameters();
};

/**
 * @brief Construct a new Parameters:: Parameters object
 *
 */
Parameters::Parameters() {}

/**
 * @brief Construct a new Parameters:: Parameters object
 *
 * @tparam Args
 * @param args
 */
template <class... Args>
Parameters::Parameters(const Args&... args) {
  this->vect_params_ = std::vector<Parameter>{args...};
}

/**
 * @brief add a new parameters
 *
 * @param param parameter to add
 */
void Parameters::add(const Parameter& param) { this->vect_params_.emplace_back(param); }

/**
 * @brief get double value of a parameter by name
 *
 * @param name name of the parameter
 * @return double double value of the parameter
 */
double Parameters::get_parameter_value(const std::string& name) const {
  auto value = this->get_val(name);
  const auto lowest_float = std::numeric_limits<float>::lowest();

  if (value > lowest_float) {
    return value;
  } else {
    throw std::runtime_error("Parameter " + name + " not found");
  }
}

/**
 * @brief get double value of a parameter by name if defined, or its default value
 *
 * @param name name of the parameter
 * @return double double value of the parameter
 */
double Parameters::get_parameter_value_or_default(const std::string& name,
                                                  const double& default_value) const {
  auto value = this->get_val(name);
  const auto lowest_float = std::numeric_limits<float>::lowest();

  if (value > lowest_float) {
    return value;
  } else {
    return default_value;
  }
}

/**
 * @brief get double value of a parameter by name
 *
 * @param name name of the parameter
 * @return double double value of the parameter
 */
double Parameters::get_val(const std::string& name) const {
  auto value = std::numeric_limits<double>::lowest();

  for (const auto& p : this->vect_params_) {
    auto pn = p.getName();
    if (pn == name) {
      value = std::get<double>(p.getValue());
    }
  }
  return value;
}

/**
 * @brief  get boolean option of a parameter by name
 *
 * @param name
 * @return true
 * @return false
 */
bool Parameters::get_option_value(const std::string& name) const {
  bool active_option = false;

  for (const auto& p : this->vect_params_) {
    auto pn = p.getName();
    if (pn == name) {
      active_option = std::get<bool>(p.getValue());
    }
  }

  return active_option;
}

/**
 * @brief transform list of parameters into a map<string,double>
 *
 * @return std::map<std::string, double>
 */
std::map<std::string, double> Parameters::getMapParameters() const {
  std::map<std::string, double> map_par;
  for (auto p : this->vect_params_) {
    auto name = p.getName();
    auto value = std::get<double>(p.getValue());
    map_par.try_emplace(name, value);
  }
  return map_par;
}

/**
 * @brief Destroy the Parameters:: Parameters object
 *
 */
Parameters::~Parameters() {}
