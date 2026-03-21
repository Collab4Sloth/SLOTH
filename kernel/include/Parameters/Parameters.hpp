/**
 * @file Parameters.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Parameter class used to build and manage calculation parameter
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
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

#pragma once

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

#include "Parameters/Parameter.hpp"

template <typename T>
concept ParameterType = std::same_as<std::remove_cvref_t<T>, Parameter>;

class Parameters {
 private:
  std::vector<Parameter> vect_params_;

 public:
  template <typename... Args>
    requires(ParameterType<Args> && ...)
  explicit Parameters(Args&&... args);
  explicit Parameters(const std::vector<Parameter>& vect_params);
  virtual ~Parameters() = default;

  std::optional<Parameter> get_parameter(const std::string& name) const;
  template <typename T>
  T get_param_value_or_default(const std::string& name, T default_value) const;

  template <typename T>
  T get_param_value(const std::string& name) const;
  std::vector<Parameter> get_vector() const;

  bool has_parameter(const std::string& param_name) const;

  int get_size() const;
  Parameters operator+(const Parameters& params) const;
  Parameters operator+(const Parameter& param) const;
  Parameters operator-(const Parameter& param) const;
};

#include "Parameters/Parameters.tpp"
