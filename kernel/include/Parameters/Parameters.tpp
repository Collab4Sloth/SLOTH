/**
 * @file Parameters.tpp
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

/**
 * @brief Construct a new Parameters:: Parameters object
 *
 * @tparam Args
 * @param args
 */
template <typename... Args>
  requires(ParameterType<Args> && ...)
Parameters::Parameters(Args&&... args) : vect_params_{std::forward<Args>(args)...} {}

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