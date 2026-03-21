/**
 * @file Parameter.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Parameter class
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
#include <string>
#include <type_traits>
#include <typeinfo>
#include <variant>

#include "Options/Options.hpp"

using param_type = std::variant<int, double, std::string, bool, MapStringDouble,
                                vTuple2StringDouble, Map2String2Double, MapString2Double, vString,
                                vInt, vDouble, vTupleStringInt, vTupleStringString>;

class Parameter {
 private:
  std::string name_;
  param_type value_;
  std::string description_;

 public:
  Parameter(const std::string& name, param_type value);

  Parameter(const std::string& name, param_type value, const std::string& description);

  std::string get_name() const;
  std::string get_description() const;

  void print() const;

  // Member function that doesn't modify the return value
  // Type is specifically mentioned despite of auto
  auto get_value() const -> param_type;

  virtual ~Parameter() = default;
};
