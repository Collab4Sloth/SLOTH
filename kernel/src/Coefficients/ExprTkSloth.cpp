/**
 * @file ExprTkCoefficient.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class definining a coefficient from an analytical formula
 * @remark Depends on ExprTk library (MIT license)
 * @version 0.1
 * @date 2025-11-06

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

#include "Coefficients/ExprTkSloth.hpp"

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
 * @brief Evaluation of the analytical expression
 *
 * @tparam T
 * @tparam Args
 * @param args
 * @return T
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
ExprTkCoefficient::F() {
  auto func = [&](const std::span<const double>& values,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    if (this->variable_names_.size() > 0) {
      if (values.size() != this->variable_names_.size()) {
        throw std::runtime_error("Number of variables not consistent with analytical expression");
      }
      for (size_t i = 0; i < values.size(); ++i) {
        this->reference_map_[this->variable_names_[i]] = values[i];
      }
    }
    return this->expression_parser_.value();
  };
  return func;
}

/**
 * @brief Function used to evaluate the gradient
 *
 * @return std::function<std::vector<double>(const std::span<const double>&)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
ExprTkCoefficient::GradientF() {
  auto func = [&](const std::span<const double>& values,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    const auto n = values.size();
    std::vector<double> gradient(n, 0.0);
    if (this->gradient_expression_parser_.size() > 0) {
      if (values.size() != this->variable_names_.size()) {
        throw std::runtime_error("Number of variables not consistent with analytical expression");
      }
      for (size_t i = 0; i < values.size(); ++i) {
        this->reference_map_[this->variable_names_[i]] = values[i];
      }
      for (size_t i = 0; i < values.size(); ++i) {
        gradient[i] = this->gradient_expression_parser_[i].value();
      }
    }
    return gradient;
  };
  return func;
}

/**
 * @brief Function used to evaluate the Hessian
 *
 * @return std::function<std::vector<double>(const std::span<const double>&)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
ExprTkCoefficient::HessianF() {
  auto func = [&](const std::span<const double>& values,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    const auto n = values.size();
    std::vector<double> hessian(n * n, 0.0);
    if (this->hessian_expression_parser_.size() > 0) {
      if (n != this->variable_names_.size()) {
        throw std::runtime_error("Number of variables not consistent with analytical expression");
      }
      for (size_t i = 0; i < n; ++i) {
        this->reference_map_[this->variable_names_[i]] = values[i];
      }
      for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
          const auto index = i * n + j;
          hessian[index] = this->hessian_expression_parser_[index].value();
        }
      }
    }
    return hessian;
  };
  return func;
}

/**
 * @brief Error manager when parsing the analytical coefficient
 *
 * @tparam T
 * @param expression Analytical expression of the coefficient
 */

void ExprTkCoefficient::ParserCoefficientError(const std::string& expression) {
  std::ostringstream oss;
  oss << "Failed to compile expression: '" << expression << "'\n";

  const std::size_t error_count = parser_.error_count();
  for (std::size_t i = 0; i < error_count; ++i) {
    auto error = parser_.get_error(i);
    oss << "Error #" << i + 1 << " : " << error.diagnostic
        << " at position: " << error.token.position << "'\n ";
  }

  throw std::runtime_error(oss.str());
}

/**
 * @brief Build the gradient
 *
 */
void ExprTkCoefficient::build_gradient() {
  for (size_t i = 0; i < gradient_math_expression_.size(); i++) {
    exprtk::expression<double> grad_parser;
    grad_parser.register_symbol_table(this->symbol_table_);
    this->gradient_expression_parser_.emplace_back(std::move(grad_parser));

    if (!this->parser_.compile(this->gradient_math_expression_[i],
                               this->gradient_expression_parser_[i])) {
      this->ParserCoefficientError(this->gradient_math_expression_[i]);
    }
  }
}

/**
 * @brief Build the Hessian
 *
 */
void ExprTkCoefficient::build_hessian() {
  for (size_t i = 0; i < hessian_math_expression_.size(); i++) {
    exprtk::expression<double> hess_parser;
    hess_parser.register_symbol_table(this->symbol_table_);
    this->hessian_expression_parser_.emplace_back(std::move(hess_parser));

    if (!this->parser_.compile(this->hessian_math_expression_[i],
                               this->hessian_expression_parser_[i])) {
      this->ParserCoefficientError(this->hessian_math_expression_[i]);
    }
  }
}

/**
 * @brief Destroy the ExprTkCoefficient::ExprTkCoefficient object
 *
 * @tparam T
 */

ExprTkCoefficient::~ExprTkCoefficient() {}
