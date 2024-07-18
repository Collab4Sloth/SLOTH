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
#include <variant>
#pragma once

using param_type = std::variant<int, double, std::string, bool>;
class Parameter {
 private:
  std::string name_;
  std::string description_;
  param_type value_;

 public:
  Parameter(const std::string& name, param_type value)
      : name_(name), value_(value), description_("") {}

  Parameter(const std::string& name, param_type value, const std::string& description)
      : name_(name), value_(value), description_(description) {}

  std::string get_name() const;
  std::string get_description() const;

  void print() const;

  // Fonction membre qui ne modifie la valeur de retour
  // Le type est spécifiquement indiqué malgré le auu
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
  std::cout << param_name << " = " << std::get<std::string>(param_value) << std::endl;
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
                      std::is_same_v<T, std::string>) {
          return arg;
        } else {
          throw std::runtime_error("Unsupported type");
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
