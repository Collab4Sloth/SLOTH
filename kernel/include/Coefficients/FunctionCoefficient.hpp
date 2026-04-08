/**
 * @file FunctionCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class defining analytical expression, first and second derivatives of a Coefficient
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
#include <any>
#include <cmath>
#include <concepts>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <utility>
#include <vector>

/**
 * @brief Base class defining analytical expression, first and second derivatives of a Coefficient
 *
 */
class FunctionCoefficient {
 protected:
  double time_;
  virtual std::function<
      double(const std::span<const double>& values, const std::span<const double>& exp_values,
             const std::span<const double>& aux_values, const unsigned int dimension)>
  F() = 0;
  virtual std::function<std::vector<double>(
      const std::span<const double>& values, const std::span<const double>& exp_values,
      const std::span<const double>& aux_values, const unsigned int dimension)>
  GradientF() = 0;
  virtual std::function<std::vector<double>(
      const std::span<const double>& values, const std::span<const double>& exp_values,
      const std::span<const double>& aux_values, const unsigned int dimension)>
  HessianF() = 0;

 public:
  FunctionCoefficient();
  virtual ~FunctionCoefficient() = default;

  // Implicit/Explicit without auxiliary variables
  double eval_f(const std::span<const double>& input_vector,
                std::optional<int> dimension = std::nullopt);
  double eval_gradient(const int i, const std::span<const double>& input_vector,
                       std::optional<int> dimension = std::nullopt);
  double eval_hessian(const int i, const int j, const std::span<const double>& input_vector,
                      std::optional<int> dimension = std::nullopt);

  // Semi-implicit without auxiliary variables
  // and Implicit/Explicit with auxiliary variables
  double eval_f(const std::span<const double>& input_vector,
                const std::span<const double>& auxiliary_vector,
                std::optional<int> dimension = std::nullopt);
  double eval_gradient(const int i, const std::span<const double>& input_vector,
                       const std::span<const double>& auxiliary_vector,
                       std::optional<int> dimension = std::nullopt);
  double eval_hessian(const int i, const int j, const std::span<const double>& input_vector,
                      const std::span<const double>& auxiliary_vector,
                      std::optional<int> dimension = std::nullopt);

  // Semi-implicit with auxiliary variables
  double eval_f(const std::span<const double>& input_vector,
                const std::span<const double>& exp_input_vector,
                const std::span<const double>& auxiliary_vector,
                std::optional<int> dimension = std::nullopt);
  double eval_gradient(const int i, const std::span<const double>& input_vector,
                       const std::span<const double>& exp_input_vector,
                       const std::span<const double>& auxiliary_vector,
                       std::optional<int> dimension = std::nullopt);
  double eval_hessian(const int i, const int j, const std::span<const double>& input_vector,
                      const std::span<const double>& exp_input_vector,
                      const std::span<const double>& auxiliary_vector,
                      std::optional<int> dimension = std::nullopt);

  virtual void set_time(double time);
};
