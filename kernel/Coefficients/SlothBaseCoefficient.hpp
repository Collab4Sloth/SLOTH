/**
 * @file SlothBaseCoefficient.hpp
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

#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "Glossary/Glossary.hpp"
#include "exprtk.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

/**
 * @brief Base class definining a coefficient from an analytical formula
 *
 * @tparam T
 */
template <class T>
class SlothBaseCoefficient {
 private:
  GlossaryQuantity coefficient_type_;
  exprtk::expression<T> expression_parser_;
  exprtk::symbol_table<T> symbol_table_;
  std::vector<std::string> variable_names_;
  std::vector<std::string> constant_names_;
  std::unordered_map<std::string, T> reference_map_;

  exprtk::parser<T> parser_;
  void ParserCoefficientError(const std::string& expression);

 public:
  SlothBaseCoefficient(GlossaryQuantity type, const std::string& expression,
                       const std::vector<std::string>& variable_names);
  SlothBaseCoefficient(GlossaryQuantity type, const std::string& expression,
                       const std::vector<std::string>& variable_names,
                       const std::vector<std::string>& constant_names);
  ~SlothBaseCoefficient();

  template <typename... Args>
  T evaluate(Args... args);

  template <typename... VarArgs, typename... ConstArgs>
  T evaluate(const std::tuple<VarArgs...>& vars, const std::tuple<ConstArgs...>& consts);

  GlossaryType get_type() const;
};

/**
 * @brief Construct a new SlothBaseCoefficient<T>::SlothBaseCoefficient object
 *
 * @tparam T
 * @param expression Analytical expression to parse
 * @param variable_names List of variables in the analytical expression (string type)
 */
template <class T>
SlothBaseCoefficient<T>::SlothBaseCoefficient(GlossaryQuantity type, const std::string& expression,
                                              const std::vector<std::string>& variable_names)
    : coefficient_type_(type), variable_names_(variable_names) {
  // Variables
  for (const auto& var : this->variable_names_) {
    this->reference_map_[var] = 0.0;  // initialize all to 0
    this->symbol_table_.add_variable(var, this->reference_map_[var]);
  }

  // Constant parameters
  this->symbol_table_.add_constants();
  this->expression_parser_.register_symbol_table(this->symbol_table_);

  if (!this->parser_.compile(expression, this->expression_parser_)) {
    this->ParserCoefficientError(expression);
  }
}

/**
 * @brief Construct a new SlothBaseCoefficient<T>::SlothBaseCoefficientobject with constants
 * @remark Constant are treated as variable for dynamic aspects
 *
 * @tparam T
 * @param expression
 * @param variable_names
 * @param constant_names
 */
template <class T>
SlothBaseCoefficient<T>::SlothBaseCoefficient(GlossaryQuantity type, const std::string& expression,
                                              const std::vector<std::string>& variable_names,
                                              const std::vector<std::string>& constant_names)
    : coefficient_type_(type), variable_names_(variable_names), constant_names_(constant_names) {
  // Variables
  for (const auto& var : this->variable_names_) {
    this->reference_map_[var] = 0.0;
    this->symbol_table_.add_variable(var, this->reference_map_[var]);
  }

  // Constant parameters
  for (const auto& var : this->constant_names_) {
    this->reference_map_[var] = 0.0;
    this->symbol_table_.add_variable(var, this->reference_map_[var]);
  }
  // this->symbol_table_.add_constants();
  this->expression_parser_.register_symbol_table(this->symbol_table_);

  if (!this->parser_.compile(expression, this->expression_parser_)) {
    this->ParserCoefficientError(expression);
  }
}

/**
 * @brief Evaluation of the analytical expression without constant parameters
 *
 * @tparam T
 * @tparam Args
 * @param args
 * @return T
 */
template <class T>
template <typename... Args>
T SlothBaseCoefficient<T>::evaluate(Args... args) {
  std::vector<T> values = {static_cast<T>(args)...};
  if (values.size() != this->variable_names_.size()) {
    throw std::runtime_error("Number of variables not consistent with analytical expression");
  }
  for (size_t i = 0; i < values.size(); ++i) {
    this->reference_map_[this->variable_names_[i]] = values[i];
  }

  return this->expression_parser_.value();
}

/**
 * @brief Evaluation of the analytical expression with constant parameters
 *
 * @tparam T
 * @tparam Args
 * @tparam Args
 * @param args
 * @param args2
 * @return T
 */
template <class T>
template <typename... VarArgs, typename... ConstArgs>
T SlothBaseCoefficient<T>::evaluate(const std::tuple<VarArgs...>& vars,
                                    const std::tuple<ConstArgs...>& consts) {
  // Check consistency of tuple sizes
  if constexpr (sizeof...(VarArgs) != 0) {
    if (sizeof...(VarArgs) != variable_names_.size()) {
      throw std::runtime_error("Number of variables not consistent with analytical expression");
    }
  }

  if constexpr (sizeof...(ConstArgs) != 0) {
    if (sizeof...(ConstArgs) != constant_names_.size()) {
      throw std::runtime_error("Number of constants not consistent with analytical expression");
    }
  }
  // Variables values
  std::apply(
      [&](auto&&... vals) {
        size_t i = 0;
        ((reference_map_[variable_names_[i++]] = vals), ...);
      },
      vars);

  // Constant values treated as variables for dynamic aspects
  std::apply(
      [&](auto&&... vals) {
        size_t i = 0;
        ((reference_map_[constant_names_[i++]] = vals), ...);
      },
      consts);

  return expression_parser_.value();
}

/**
 * @brief Error manager when parsing the analytical coefficient
 *
 * @tparam T
 * @param expression Analytical expression of the coefficient
 */
template <class T>
void SlothBaseCoefficient<T>::ParserCoefficientError(const std::string& expression) {
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
 * @brief Return the GlobalQuantity associated with the coefficient
 *
 * @tparam T
 * @return GlossaryQuantity
 */
template <class T>
GlossaryType SlothBaseCoefficient<T>::get_type() const {
  return this->coefficient_type_.type;
}

/**
 * @brief Destroy the SlothBaseCoefficient<T>::SlothBaseCoefficient object
 *
 * @tparam T
 */
template <class T>
SlothBaseCoefficient<T>::~SlothBaseCoefficient() {}