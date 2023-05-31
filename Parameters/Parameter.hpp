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
#include <set>
#include <string>
#include <variant>
#pragma once

using var = std::variant<int, double, std::string, bool>;
class Parameter {
 private:
  std::string name;
  var value;

 public:
  Parameter(std::string name, var value) : name(name), value(value) {
    // TODO(ci230846) : définir une méthode de check sur base std::visit et overload
  }

  /** Method used to get the name of the parameter
   *  return name of the parameter of type string
   */
  std::string getName() const { return name; }  // end of getName

  /** Method used to get the value of the parameter
   *  return value of the parameter of any type (see variant)
   */
  var getValue() const { return value; }  // end of getValue

  void pprint() {
    std::visit([](const auto& x) { std::cout << x << std::endl; }, value);
  }

  ~Parameter() {}
};
