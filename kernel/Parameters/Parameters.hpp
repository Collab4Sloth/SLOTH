/*
 * Copyright © CEA 2022
 *
 * \brief Parameter class used to build and manage calculation parameter
 *
 * \file Parameter.hpp
 * \author ci230846
 * \date 19/01/2022
 */
#include <algorithm>
#include <any>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "Parameters/Parameter.hpp"
#pragma once

class Parameters {
 private:
  std::vector<Parameter> vect_params_;

 public:
  template <typename... Args>
  explicit Parameters(Args&&... args);
  explicit Parameters(const std::vector<Parameter>& vect_params);

  std::optional<Parameter> get_parameter(const std::string& name) const;
  template <typename T>
  T get_param_value_or_default(const std::string& name, T default_value) const;

  template <typename T>
  T get_param_value(const std::string& name) const;
  std::vector<Parameter> get_vector() const;

  bool has_parameter(const std::string& param_name) const;

  int get_size() const;

  /**
   * @brief Add two Parameters with a priority given to the second Parameters to overwrite
   * Parameter objets in the first Parameters if already exist.
   *
   * @param params
   * @return Parameters
   */
  Parameters operator+(const Parameters& params) const {
    auto initial_vect = this->get_vector();
    auto add_vect = params.get_vector();
    std::vector<Parameter> sum;
    if (!std::empty(initial_vect)) {
      sum.insert(sum.begin(), initial_vect.begin(), initial_vect.end());
    }
    if (!std::empty(add_vect)) {
      for (const auto& add_param : add_vect) {
        sum.erase(std::remove_if(
                      sum.begin(), sum.end(),
                      [&add_param](const auto& v) { return v.get_name() == add_param.get_name(); }),
                  sum.end());
        sum.emplace_back(add_param);
      }
    }
    return Parameters(sum);
  }
  /**
   * @brief Add a parameter inside a Parameters. If param already exists, the given parameters
   * overwrites the existing one.
   *
   * @param param
   * @return Parameters
   */
  Parameters operator+(const Parameter& param) const {
    auto initial_vect = this->get_vector();
    // Remove a parameter inside a Parameters if exists
    initial_vect.erase(
        std::remove_if(initial_vect.begin(), initial_vect.end(),
                       [&param](const auto& v) { return v.get_name() == param.get_name(); }),
        initial_vect.end());

    initial_vect.emplace_back(param);

    return Parameters(initial_vect);
  }
  /**
   * @brief Remove a parameter inside a Parameters if exists
   *
   * @param param
   * @return Parameters
   */
  Parameters operator-(const Parameter& param) const {
    auto initial_vect = this->get_vector();
    initial_vect.erase(
        std::remove_if(initial_vect.begin(), initial_vect.end(),
                       [&param](const auto& v) { return v.get_name() == param.get_name(); }),
        initial_vect.end());

    return Parameters(initial_vect);
  }

  ~Parameters() = default;
};

/**
 * @brief Construct a new Parameters:: Parameters object
 *
 * @tparam Args
 * @param args
 */
template <typename... Args>
Parameters::Parameters(Args&&... args) {
  this->vect_params_ = std::vector<Parameter>{std::forward<Args>(args)...};
}

/**
 * @brief Construct a new Parameters:: Parameters object
 *
 * @param vect_params
 */
DEBILE_INLINE Parameters::Parameters(const std::vector<Parameter>& vect_params) : vect_params_(vect_params) {}

/**
 * @brief Return a vector of Parameter
 *
 * @tparam Args
 * @return std::vector<Parameter>
 */
DEBILE_INLINE std::vector<Parameter> Parameters::get_vector() const { return this->vect_params_; }

/**
 * @brief Search a parameter by its name
 *
 * @param param_name
 * @return true
 * @return false
 */
DEBILE_INLINE bool Parameters::has_parameter(const std::string& param_name) const {
  return std::find_if(vect_params_.begin(), vect_params_.end(), [&](const auto& p) {
           return p.get_name() == param_name;
         }) != vect_params_.end();
}

/**
 * @brief Return the size of a Parameters
 *
 * @tparam Args
 * @return in
 */
DEBILE_INLINE int Parameters::get_size() const { return this->vect_params_.size(); }

/**
 * @brief Return the parameter associated to the given name
 *
 * @param name
 * @return Parameter
 */
DEBILE_INLINE std::optional<Parameter> Parameters::get_parameter(const std::string& name) const {
  for (const auto& param : this->vect_params_) {
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
      // int rank = mfem::Mpi::WorldRank();
      // if (rank == 0) {
      //   std::string warn_mess =
      //       "##############\n Parameter named " + name + "  overloaded.\n##############";
      //   mfem::mfem_warning(warn_mess.c_str());
      // }

      return std::get<T>(value);
    } catch (const std::bad_variant_access&) {
      std::string error_mess = "##############\n Parameter named " + name +
                               ": Invalid conversion. Please check the data.\n##############";
      mfem::mfem_error(error_mess.c_str());
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
      std::string error_mess = "##############\n Parameter named " + name +
                               " Invalid conversion. Please check the data.\n##############";
      mfem::mfem_error(error_mess.c_str());
    }
  } else {
    std::string error_mess = "##############\n Parameter named " + name +
                             " Not Found. Please check the data.\n##############";
    mfem::mfem_error(error_mess.c_str());
  }
}
