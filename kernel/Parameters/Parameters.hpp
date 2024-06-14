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
#include <optional>

#include "Parameters/Parameter.hpp"
#pragma once

class Parameters {
 private:
  std::vector<Parameter> vect_params_;

 public:
  template <class... Args>
  Parameters(Args&&... args);
  std::optional<Parameter> get_parameter(const std::string& name) const;
  template <typename T>
  T get_param_value_or_default(const std::string& name, T default_value) const;
  template <typename T>
  T get_param_value(const std::string& name) const;

  ~Parameters();
};

/**
 * @brief Construct a new Parameters:: Parameters object
 *
 * @tparam Args
 * @param args
 */
template <class... Args>
Parameters::Parameters(Args&&... args) {
  this->vect_params_ = std::vector<Parameter>{std::forward<Args>(args)...};
}

/**
 * @brief Return the parameter associated to the given name
 *
 * @param name
 * @return Parameter
 */
std::optional<Parameter> Parameters::get_parameter(const std::string& name) const {
  for (const auto& param : vect_params_) {
    if (param.get_name() == name) {
      return param;
    }
  }
  return std::nullopt;
}
/**
 * @brief Return the value of the parameter associated to the given name or, use th edefault given,
 * value
 * @param name
 * @param default_value
 * @return std::variant<int, double, std::string>
 */
template <typename T>
T Parameters::get_param_value_or_default(const std::string& name, T default_value) const {
  auto xx = this->get_parameter(name);
  if (xx) {
    const Parameter& pp = xx.value();
    try {
      auto value = pp.get_value();
      return std::get<T>(value);
    } catch (const std::bad_variant_access&) {
      std::string error_mess = name + "Invalid conversion. Please check the data.";
      throw std::runtime_error(error_mess);
    }
  } else {
    return default_value;
  }
}

/**
 * @brief Return the value of the parameter associated to the given name
 *
 * @param name
 * @return std::variant<int, double, std::string>
 */
template <typename T>
T Parameters::get_param_value(const std::string& name) const {
  auto xx = this->get_parameter(name);
  if (xx) {
    const Parameter& pp = xx.value();
    try {
      auto value = pp.get_value();
      return std::get<T>(value);
    } catch (const std::bad_variant_access&) {
      std::string error_mess = name + " Invalid conversion. Please check the data.";
      throw std::runtime_error(error_mess);
    }
  } else {
    std::string error_mess = name + "Not Found. Please check the data.";
    throw std::runtime_error(error_mess);
  }
}

/**
 * @brief Destroy the Parameters:: Parameters object
 *
 */
Parameters::~Parameters() {}
