/**
 * @file Coefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class definining a coefficient from an analytical formula
 * @version 0.1
 * @date 2025-11-06
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

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "Coefficients/FileCoefficient.hpp"
#include "Coefficients/FunctionCoefficient.hpp"
#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Glossary/Glossary.hpp"
#include "Options/TimeOptions.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Base class definining a coefficient from an analytical formula
 *
 */
class Coefficient : public SlothBaseCoefficient {
 public:
  /**
   * @brief Builds a coefficient from a non-file FunctionCoefficient passed as an rvalue.
   *        The object is moved into a shared_ptr<FunctionCoefficient>.
   *
   * @tparam T Concrete FunctionCoefficient type (FileCoefficient excluded)
   * @param type Quantity associated with this coefficient
   * @param scheme Discretization scheme used by this coefficient
   * @param coef Function coefficient to wrap
   */
  template <class T>
    requires std::derived_from<std::remove_cvref_t<T>, FunctionCoefficient> &&
             (!std::derived_from<std::remove_cvref_t<T>, FileCoefficient>)
  Coefficient(GlossaryQuantity type, Scheme scheme, T&& coef)
      : Coefficient(type, scheme, std::make_shared<std::remove_cvref_t<T>>(std::forward<T>(coef))) {
  }

  /**
   * @brief Builds a coefficient from a non-file FunctionCoefficient passed as a const lvalue.
   *        The object is copied into a shared_ptr<FunctionCoefficient>.
   *
   * @tparam T Concrete FunctionCoefficient type (FileCoefficient excluded)
   * @param type Quantity associated with this coefficient
   * @param scheme Discretization scheme used by this coefficient
   * @param coef Function coefficient to wrap
   */
  template <class T>
    requires std::derived_from<T, FunctionCoefficient> && (!std::derived_from<T, FileCoefficient>)
  Coefficient(GlossaryQuantity type, Scheme scheme, const T& coef)
      : Coefficient(type, scheme, std::make_shared<T>(coef)) {}

  /**
   * @brief Builds a coefficient from a FileCoefficient (scheme is implicitly Scheme::File).
   *
   * @tparam T Concrete FileCoefficient type
   * @param type Quantity associated with this coefficient
   * @param coef File coefficient to wrap (moved if rvalue, copied otherwise)
   */
  template <class T>
    requires std::derived_from<std::remove_cvref_t<T>, FileCoefficient>
  Coefficient(GlossaryQuantity type, T&& coef)
      : Coefficient(type, std::make_shared<std::remove_cvref_t<T>>(std::forward<T>(coef))) {}

  /**
   * @brief Builds a coefficient from a FileCoefficient with Glossary::Default as quantity
   *        (scheme is implicitly Scheme::File).
   *
   * @tparam T Concrete FileCoefficient type
   * @param coef File coefficient to wrap (moved if rvalue, copied otherwise)
   */
  template <class T>
    requires std::derived_from<std::remove_cvref_t<T>, FileCoefficient>
  explicit Coefficient(T&& coef)
      : Coefficient(std::make_shared<std::remove_cvref_t<T>>(std::forward<T>(coef))) {}

  /**
   * @brief Builds a coefficient from a shared FunctionCoefficient.
   *
   * @param type Quantity associated with this coefficient
   * @param scheme Discretization scheme used by this coefficient
   * @param coef Shared function coefficient
   */
  Coefficient(GlossaryQuantity type, Scheme scheme, std::shared_ptr<FunctionCoefficient> coef);

  /**
   * @brief Builds a coefficient from a shared FileCoefficient (scheme is Scheme::File).
   *
   * @param type Quantity associated with this coefficient
   * @param coef Shared file coefficient
   */
  Coefficient(GlossaryQuantity type, std::shared_ptr<FileCoefficient> coef);

  /**
   * @brief Builds a coefficient from a shared FileCoefficient with Glossary::Default as quantity
   *        (scheme is Scheme::File).
   *
   * @param coef Shared file coefficient
   */
  explicit Coefficient(std::shared_ptr<FileCoefficient> coef);

  /**
   * @brief Builds a constant coefficient.
   *
   * @param type Quantity associated with this coefficient
   * @param coef Constant value
   */
  Coefficient(GlossaryQuantity type, double coef);

  virtual ~Coefficient() = default;
};
