/**
 * @file ExprTkCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class definining a coefficient from an analytical formula
 * @remark Depends on ExprTk library (MIT license)
 * @version 0.1
 * @date 2025-11-06
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
#pragma once

#include <iostream>
#include <span>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Coefficients/FunctionCoefficient.hpp"
#include "Glossary/Glossary.hpp"
#include "exprtk.hpp"  // NOLINT [no include the directory when naming exprtk include file]
#include "mfem.hpp"    // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Base class definining a coefficient from an analytical formula
 *
 * @tparam T
 */

class ExprTkCoefficient : public FunctionCoefficient {
 private:
  exprtk::parser<double> parser_;
  exprtk::symbol_table<double> symbol_table_;
  exprtk::expression<double> expression_parser_;
  std::vector<exprtk::expression<double>> gradient_expression_parser_;
  std::vector<exprtk::expression<double>> hessian_expression_parser_;

  std::vector<std::string> variable_names_;
  std::string math_expression_;

  std::vector<std::string> gradient_math_expression_;
  std::vector<std::string> hessian_math_expression_;

  std::unordered_map<std::string, double> reference_map_;

  void ParserCoefficientError(const std::string& expression);

  template <class... Args>
  void build_functional(Args&&... function_variables_names);
  void build_gradient();
  void build_hessian();

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const unsigned int dimension)>
  F() override final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  GradientF() override final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  HessianF() override final;

 public:
  template <class... Args>
    requires((sizeof...(Args) > 0) && (std::convertible_to<Args, std::string_view> && ...))
  explicit ExprTkCoefficient(Args&&... function_variable_names);

  template <class... Args>
    requires((sizeof...(Args) > 0) && (std::convertible_to<Args, std::string_view> && ...))
  explicit ExprTkCoefficient(std::vector<std::string> grad_functions,
                             Args&&... function_variable_names);

  template <class... Args>
    requires((sizeof...(Args) > 0) && (std::convertible_to<Args, std::string_view> && ...))
  ExprTkCoefficient(std::vector<std::string> hess_functions,
                    std::vector<std::string> grad_functions, Args&&... function_variable_names);
  ~ExprTkCoefficient();

  GlossaryType get_type() const;
};

/**
 * @brief Construct a new ExprTkCoefficient::ExprTkCoefficient object
 *
 * @tparam T
 * @param expression Analytical expression to parse
 * @param variable_names List of variables in the analytical expression (string type)
 */
template <class... Args>
  requires((sizeof...(Args) > 0) && (std::convertible_to<Args, std::string_view> && ...))
ExprTkCoefficient::ExprTkCoefficient(Args&&... function_variable_names) : FunctionCoefficient() {
  MFEM_VERIFY(sizeof...(function_variable_names) != 0,
              "Error while defining Coefficient. Please check your data");

  ///////////////////
  // Functional
  ///////////////////
  this->build_functional(function_variable_names...);
}

template <class... Args>
  requires((sizeof...(Args) > 0) && (std::convertible_to<Args, std::string_view> && ...))
ExprTkCoefficient::ExprTkCoefficient(std::vector<std::string> grad_functions,
                                     Args&&... function_variable_names)
    : FunctionCoefficient(), gradient_math_expression_(grad_functions) {
  MFEM_VERIFY(sizeof...(function_variable_names) != 0,
              "Error while defining Coefficient. Please check your data");

  ///////////////////
  // Functional
  ///////////////////
  this->build_functional(function_variable_names...);

  ///////////////////
  // Gradient
  ///////////////////
  this->build_gradient();
}

template <class... Args>
  requires((sizeof...(Args) > 0) && (std::convertible_to<Args, std::string_view> && ...))
ExprTkCoefficient::ExprTkCoefficient(std::vector<std::string> hess_functions,
                                     std::vector<std::string> grad_functions,
                                     Args&&... function_variable_names)
    : FunctionCoefficient(),
      gradient_math_expression_(grad_functions),
      hessian_math_expression_(hess_functions) {
  MFEM_VERIFY(sizeof...(function_variable_names) != 0,
              "Error while defining Coefficient. Please check your data");

  ///////////////////
  // Functional
  ///////////////////

  this->build_functional(function_variable_names...);

  ///////////////////
  // Gradient
  ///////////////////
  this->build_gradient();

  ///////////////////
  // Hessian
  ///////////////////
  this->build_hessian();
}

/**
 * @brief Build the function
 *
 * @tparam Args
 * @param function_variables_names
 */
template <class... Args>
void ExprTkCoefficient::build_functional(Args&&... function_variable_names) {
  this->math_expression_ = "";
  // Constant mathematical expression
  if constexpr (sizeof...(function_variable_names) == 1) {
    this->math_expression_ = {std::forward<Args>(function_variable_names)...};
    this->variable_names_.clear();
  } else {
    // Variable mathematical expression
    this->variable_names_ = {std::forward<Args>(function_variable_names)...};

    this->math_expression_ = this->variable_names_[0];
    if (!this->variable_names_.empty()) {
      this->variable_names_.erase(this->variable_names_.begin());
    }
  }
  // Variables
  for (const auto& var : this->variable_names_) {
    this->reference_map_[var] = 0.0;  // initialize all to 0
    this->symbol_table_.add_variable(var, this->reference_map_[var]);
  }

  // Constant parameters
  this->symbol_table_.add_constants();
  this->expression_parser_.register_symbol_table(this->symbol_table_);

  if (!this->parser_.compile(this->math_expression_, this->expression_parser_)) {
    this->ParserCoefficientError(this->math_expression_);
  }
}
