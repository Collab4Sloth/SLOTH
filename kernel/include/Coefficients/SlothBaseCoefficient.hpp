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
#pragma once
#include <iostream>
#include <memory>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include "Coefficients/ConstantCoefficient.hpp"
#include "Coefficients/FunctionCoefficient.hpp"
#include "Glossary/Glossary.hpp"
#include "Options/TimeOptions.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Base class defining a coefficient from an analytical formula.
 */
class SlothBaseCoefficient {
 private:
  GlossaryQuantity coefficient_quantity_;
  std::shared_ptr<FunctionCoefficient> coefficient_;
  std::vector<int> bdr_index_;

  Scheme scheme_;

 public:
  SlothBaseCoefficient(GlossaryQuantity type, Scheme scheme,
                       std::shared_ptr<FunctionCoefficient> coef);
  SlothBaseCoefficient(GlossaryQuantity type, double coef);

  virtual ~SlothBaseCoefficient() = default;

  // Scalar
  double compute();

  // Implicit/Explicit without auxiliary variables
  double compute(const std::span<const double>& values,
                 std::optional<int> dimension = std::nullopt);

  // Semi-implicit without auxiliary variables
  // and Implicit/Explicit with auxiliary variables
  double compute(const std::span<const double>& values,
                 const std::span<const double>& auxiliary_values,
                 std::optional<int> dimension = std::nullopt);

  // Semi-implicit with auxiliary variables
  double compute(const std::span<const double>& values, const std::span<const double>& exp_values,
                 const std::span<const double>& auxiliary_values,
                 std::optional<int> dimension = std::nullopt);

  // Implicit/Explicit without auxiliary variables
  double compute_gradient(const int i, const std::span<const double>& values,
                          std::optional<int> dimension = std::nullopt);

  // Semi-implicit without auxiliary variables
  // and Implicit/Explicit with auxiliary variables
  double compute_gradient(const int i, const std::span<const double>& values,
                          const std::span<const double>& auxiliary_values,
                          std::optional<int> dimension = std::nullopt);

  // Semi-implicit with auxiliary variables
  double compute_gradient(const int i, const std::span<const double>& values,
                          const std::span<const double>& exp_values,
                          const std::span<const double>& auxiliary_values,
                          std::optional<int> dimension = std::nullopt);

  // Implicit/Explicit without auxiliary variables
  double compute_hessian(const int i, const int j, const std::span<const double>& values,
                         std::optional<int> dimension = std::nullopt);

  // Semi-implicit without auxiliary variables
  // and Implicit/Explicit with auxiliary variables
  double compute_hessian(const int i, const int j, const std::span<const double>& values,
                         const std::span<const double>& auxiliary_values,
                         std::optional<int> dimension = std::nullopt);

  // Semi-implicit with auxiliary variables
  double compute_hessian(const int i, const int j, const std::span<const double>& values,
                         const std::span<const double>& exp_values,
                         const std::span<const double>& auxiliary_values,
                         std::optional<int> dimension = std::nullopt);

  GlossaryType get_type() const;
  const GlossaryQuantity get_quantity();
  unsigned int get_id() const;
  bool is_implicit() const;
  bool is_explicit() const;
  bool is_semi_implicit() const;
  bool is_scalar() const;

  void set_time(double time);
  void set_bdr_index_coef(std::vector<int> ids);
  std::vector<int> get_bdr_index_coef();
};
