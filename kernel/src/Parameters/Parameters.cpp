/**
 * @file Parameters.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Parameter class used to build and manage calculation parameter
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
#include "Parameters/Parameters.hpp"

#include <algorithm>
#include <any>
#include <concepts>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

/**
 * @brief Construct a new Parameters:: Parameters object
 *
 * @param vect_params
 */
Parameters::Parameters(const std::vector<Parameter>& vect_params) : vect_params_(vect_params) {}

/**
 * @brief Return a vector of Parameter
 *
 * @tparam Args
 * @return std::vector<Parameter>
 */
std::vector<Parameter> Parameters::get_vector() const { return this->vect_params_; }

/**
 * @brief Search a parameter by its name
 *
 * @param param_name
 * @return true
 * @return false
 */
bool Parameters::has_parameter(const std::string& param_name) const {
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
int Parameters::get_size() const { return this->vect_params_.size(); }

/**
 * @brief Return the parameter associated to the given name
 *
 * @param name
 * @return Parameter
 */
std::optional<Parameter> Parameters::get_parameter(const std::string& name) const {
  for (const auto& param : this->vect_params_) {
    if (param.get_name() == name) {
      return param;
    }
  }
  return std::nullopt;
}

/**
 * @brief Remove a parameter inside a Parameters if exists
 *
 * @param param
 * @return Parameters
 */
Parameters Parameters::operator-(const Parameter& param) const {
  auto initial_vect = this->get_vector();
  initial_vect.erase(
      std::remove_if(initial_vect.begin(), initial_vect.end(),
                     [&param](const auto& v) { return v.get_name() == param.get_name(); }),
      initial_vect.end());

  return Parameters(initial_vect);
}

/**
 * @brief Add a parameter inside a Parameters. If param already exists, the given parameters
 * overwrites the existing one.
 *
 * @param param
 * @return Parameters
 */
Parameters Parameters::operator+(const Parameter& param) const {
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
 * @brief Add two Parameters with a priority given to the second Parameters to overwrite
 * Parameter objets in the first Parameters if already exist.
 *
 * @param params
 * @return Parameters
 */
Parameters Parameters::operator+(const Parameters& params) const {
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
