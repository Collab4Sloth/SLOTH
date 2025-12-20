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

#include "FunctionCoefficient.hpp"
#include "Glossary/Glossary.hpp"
#include "Options/TimeOptions.hpp"
#include "SlothBaseCoefficient.hpp"

#ifdef SLOTH_USE_EXPRTK
#include "ExprTkSloth.hpp"
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

  virtual ~Coefficient() = default;
};
/**
 * @brief Construct a new Coefficient::Coefficient object
 *
 * @tparam Args
 * @param type
 * @param coef
 */
inline Coefficient::Coefficient(GlossaryQuantity type, Scheme scheme,
                                std::shared_ptr<FunctionCoefficient> coef)
    : SlothBaseCoefficient(type, scheme, std::move(coef)) {}

inline Coefficient::Coefficient(GlossaryQuantity type, double coef)
    : SlothBaseCoefficient(type, std::move(coef)) {}

#ifdef SLOTH_USE_EXPRTK
/**
 * @brief Construct a new Coefficient::Coefficient object from strings
 * @remark The first string is mandatory and corresponds to the function. If not constant, th other
 * strings corresponds to variables
 *
 * @tparam Args
 */
template <class... Args>
  requires((sizeof...(Args) > 0) && (std::convertible_to<Args, std::string_view> && ...))
inline Coefficient::Coefficient(GlossaryQuantity type, Scheme scheme,
                                Args&&... function_variable_names)
    : SlothBaseCoefficient(
          type, scheme,
          std::make_shared<ExprTkCoefficient>(std::forward<Args>(function_variable_names)...)) {}

/**
 * @brief Construct a new Coefficient::Coefficient object with its gradient and hessian
 * @remark The first vector is the hessian, the second the gradient. The first string corresponds to
 * the function. It is mandatory. If not constant, the other strings corresponds to variables
 *
 * @tparam Args
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
 * @brief Construct a new Coefficient::Coefficient object with its gradient
 * @remark The first vector is the gradient. The first string corresponds to the function. It is
 * mandatory. If not constant, the other strings corresponds to variables
 *
 * @tparam Args
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
