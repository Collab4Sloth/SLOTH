/**
 * @file Parameter.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-04-30
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <any>
#include <limits>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <variant>
#include <vector>

#include "Options/Options.hpp"

#pragma once

using param_type =
    std::variant<int, double, std::string, bool, MapStringDouble, vTuple2StringDouble,
                 Map2String2Double, vString, vTupleStringInt, vTupleStringString,
                 std::function<double(double)>, std::vector<std::size_t>>;
class Parameter {
 private:
  std::string name_;
  param_type value_;
  std::string description_;

 public:
  Parameter(const std::string& name, param_type value)
      : name_(name), value_(value), description_("") {}

  Parameter(const std::string& name, param_type value, const std::string& description)
      : name_(name), value_(value), description_(description) {}

  std::string get_name() const;
  std::string get_description() const;

  void print() const;

  // Member function that doesn't modify the return value
  // Type is specifically mentioned despite of auto
  auto get_value() const -> param_type;

  ~Parameter() {}
};

/**
 * @brief Print the value of the parameter
 *
 */
void Parameter::print() const {
  const auto& param_value = this->get_value();
  const auto& param_name = this->get_name();

  std::visit(
      [param_name](auto&& arg) {
        using T = std::decay_t<decltype(arg)>;

        if constexpr (std::is_same_v<T, int>) {
          SlothInfo::print(param_name, " = ", arg);
        } else if constexpr (std::is_same_v<T, double>) {
          SlothInfo::print(param_name, " = ", arg);
        } else if constexpr (std::is_same_v<T, std::string>) {
          SlothInfo::print(param_name, " = ", arg);
        } else if constexpr (std::is_same_v<T, bool>) {
          SlothInfo::print(param_name, " = ", arg);
        } else if constexpr (std::is_same_v<T, vString>) {
          for (const auto& v : arg) {
            SlothInfo::print(param_name, " = ", v);
          }
        } else if constexpr (std::is_same_v<T, vTupleStringString>) {
          for (const auto& v : arg) {
            const auto& [s1, s2] = v;
            SlothInfo::print(param_name, " = ", s1, ", ", s2);
          }
        } else if constexpr (std::is_same_v<T, vTupleStringInt>) {
          for (const auto& v : arg) {
            const auto& [s1, i2] = v;
            SlothInfo::print(param_name, " = ", s1, ", ", i2);
          }

        } else if constexpr (std::is_same_v<T, MapStringDouble>) {
          for (const auto& [key, val] : arg) {
            SlothInfo::print(param_name, " = ", key, ", ", val);
          }
        } else if constexpr (std::is_same_v<T, vTuple2StringDouble>) {
          for (const auto& [k1, k2, k3] : arg) {
            SlothInfo::print(param_name, " = [", k1, ", ", k2, ",", k3, "]");
          }
        } else if constexpr (std::is_same_v<T, Map2String2Double>) {
          for (const auto& [key, val] : arg) {
            const auto& [k1, k2] = key;
            const auto& [v1, v2] = val;
            SlothInfo::print(param_name, " = [", k1, ", ", k2, "]  [", v1, ", ", v2, "]");
          }
        } else if constexpr (std::is_same_v<T, std::function<double(double)>>) {
          SlothInfo::print(param_name,
                           " This parameter is a function, check source for more informations");
        } else if constexpr (std::is_same_v<T, std::vector<std::size_t>>) {
          for (size_t i = 0; i < arg.size(); i++) {
            SlothInfo::print(param_name, " [", i, "] = ", arg.at(i), ",");
          }
        } else {
          mfem::mfem_error("Unsupported type");
        }
      },
      this->value_);
}

/**
 * @brief Return the value of the paramter
 *
 * @return std::variant<int, double, std::string>
 */
auto Parameter::get_value() const -> param_type {
  return std::visit(
      [](auto&& arg) -> param_type {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, int> || std::is_same_v<T, double> ||
                      std::is_same_v<T, std::string> || std::is_same_v<T, bool> ||
                      std::is_same_v<T, vString> || std::is_same_v<T, vTupleStringString> ||
                      std::is_same_v<T, vTupleStringInt> || std::is_same_v<T, MapStringDouble> ||
<<<<<<< HEAD
                      std::is_same_v<T, Map2String2Double> ||
                      std::is_same_v<T, std::function<double(double)>> ||
                      std::is_same_v<T, std::vector<std::size_t>>) {
=======
                      std::is_same_v<T, vTuple2StringDouble> ||
                      std::is_same_v<T, Map2String2Double>) {
>>>>>>> master
          return arg;
        } else {
          mfem::mfem_error("Unsupported type");
        }
      },
      this->value_);
}

/**
 * @brief Return the name of the parameter
 *
 * @return std::string
 */
std::string Parameter::get_name() const { return this->name_; }

/**
 * @brief Return the description  of the parameter
 *
 * @return std::string
 */
std::string Parameter::get_description() const { return this->description_; }
