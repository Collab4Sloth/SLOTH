/**
 * @file Coefficient.hpp
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
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "Coefficients/FunctionCoefficient.hpp"
#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Glossary/Glossary.hpp"
#include "Options/TimeOptions.hpp"

#ifdef SLOTH_USE_EXPRTK
#include "ExprTkSloth.hpp"  // NOLINT [no include the directory when naming exprtk include file]
#endif
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once

/**
 * @brief Base class definining a coefficient from an analytical formula
 *
 */
class Coefficient : public SlothBaseCoefficient {
 public:
  /**
   * @brief Perfect-forwarding constructor that converts any FunctionCoefficient object into a
   * shared_ptr<FunctionCoefficient> before calling a second constructor.
   *
   * @tparam T
   * @param type
   * @param coef
   * @return requires
   */
  template <class T>
    requires std::derived_from<T, FunctionCoefficient>
  Coefficient(GlossaryQuantity type, Scheme scheme, T&& coef)
      : Coefficient(type, scheme, std::make_shared<std::decay_t<T>>(std::forward<T>(coef))) {}

  Coefficient(GlossaryQuantity type, Scheme scheme, std::shared_ptr<FunctionCoefficient> coef);
  Coefficient(GlossaryQuantity type, double coef);
  virtual ~Coefficient() = default;

#ifdef SLOTH_USE_EXPRTK
  template <class... Args>
    requires((sizeof...(Args) > 0) && (std::convertible_to<Args, std::string_view> && ...))
  explicit Coefficient(GlossaryQuantity type, Scheme scheme, Args&&... function_variable_names);

  template <class... Args>
    requires((sizeof...(Args) > 0) && (std::convertible_to<Args, std::string_view> && ...))
  Coefficient(GlossaryQuantity type, Scheme scheme, std::vector<std::string> grad_functions,
              Args&&... function_variable_names);

  template <class... Args>
    requires((sizeof...(Args) > 0) && (std::convertible_to<Args, std::string_view> && ...))
  Coefficient(GlossaryQuantity type, Scheme scheme, std::vector<std::string> hess_functions,
              std::vector<std::string> grad_functions, Args&&... function_variable_names);
#endif
};

/**
 * @brief Constructs a Coefficient object from a user-provided function coefficient.
 *
 * This constructor allows creating a coefficient with a custom FunctionCoefficient
 * implementation.
 *
 * @param type Quantity associated with this coefficient.
 * @param scheme Discretization scheme used by this coefficient.
 * @param coef Shared pointer to the function-based coefficient.
 */
inline Coefficient::Coefficient(GlossaryQuantity type, Scheme scheme,
                                std::shared_ptr<FunctionCoefficient> coef)
    : SlothBaseCoefficient(type, scheme, std::move(coef)) {}

/**
 * @brief Constructs a Coefficient object with a constant value.
 *
 * @param type Quantity associated with this coefficient.
 * @param coef Constant coefficient value.
 */
inline Coefficient::Coefficient(GlossaryQuantity type, double coef)
    : SlothBaseCoefficient(type, std::move(coef)) {}

#ifdef SLOTH_USE_EXPRTK
/**
 * @brief Constructs a Coefficient object from function and variable names.
 *
 * This constructor allows creating a coefficient directly from string expressions.
 * - The first string corresponds to the main function and is mandatory.
 * - If the function is not constant, additional strings correspond to the variable names.
 *
 * @tparam Args Variadic template parameters representing variable names (must be convertible to
 * `std::string_view`).
 * @param type Quantity associated with this coefficient.
 * @param scheme Discretization scheme used by this coefficient.
 * @param function_variable_names Names of the function and its variables.
 */
template <class... Args>
  requires((sizeof...(Args) > 0) && (std::convertible_to<Args, std::string_view> && ...))
inline Coefficient::Coefficient(GlossaryQuantity type, Scheme scheme,
                                Args&&... function_variable_names)
    : SlothBaseCoefficient(
          type, scheme,
          std::make_shared<ExprTkCoefficient>(std::forward<Args>(function_variable_names)...)) {}

/**
 * @brief Constructs a Coefficient object from expressions for both the Hessian and its gradient.
 *
 * This constructor allows specifying Hessian and gradient functions explicitly.
 * - hess_functions contains the expressions for the Hessian components.
 * - grad_functions contains the gradient functions.
 * - The first string in hess_functions or grad_functions corresponds to the main function and
 * is mandatory.
 * - Additional strings correspond to variable names used in the expressions.
 *
 * @tparam Args Variadic template parameters for additional variable names (must be convertible to
 * std::string_view).
 * @param type Quantity associated with this coefficient.
 * @param scheme Discretization scheme used by this coefficient.
 * @param hess_functions Vector of strings representing Hessian expressions.
 * @param grad_functions Vector of strings representing gradient expressions.
 * @param function_variable_names Names of the variables used in the expressions.
 */
template <class... Args>
  requires((sizeof...(Args) > 0) && (std::convertible_to<Args, std::string_view> && ...))
inline Coefficient::Coefficient(GlossaryQuantity type, Scheme scheme,
                                std::vector<std::string> hess_functions,
                                std::vector<std::string> grad_functions,
                                Args&&... function_variable_names)
    : SlothBaseCoefficient(
          type, scheme,
          std::make_shared<ExprTkCoefficient>(hess_functions, grad_functions,
                                              std::forward<Args>(function_variable_names)...)) {}

/**
 * @brief Constructs a Coefficient object from an analytical expression with its gradient.
 *
 * This constructor allows specifying the gradient functions explicitly.
 * The first element of grad_functions is the gradient itself. The first string
 * corresponds to the main function and is mandatory. If the function is not constant,
 * the remaining strings correspond to the variable names used in the expression.
 *
 * @tparam Args Variadic template parameters for additional variable names (must be convertible to
 * std::string_view).
 * @param type Quantity associated with this coefficient.
 * @param scheme Discretization scheme used by this coefficient.
 * @param grad_functions Vector of strings representing the gradient functions.
 * @param function_variable_names Names of the variables used in the expression.
 */
template <class... Args>
  requires((sizeof...(Args) > 0) && (std::convertible_to<Args, std::string_view> && ...))
inline Coefficient::Coefficient(GlossaryQuantity type, Scheme scheme,
                                std::vector<std::string> grad_functions,
                                Args&&... function_variable_names)
    : SlothBaseCoefficient(type, scheme,
                           std::make_shared<ExprTkCoefficient>(
                               grad_functions, std::forward<Args>(function_variable_names)...)) {}
#endif
