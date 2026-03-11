/**
 *
 * Copyright CEA (C) 2026
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

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <span>
#include <vector>

#include "Options/PhysicalPropertiesOptions.hpp"

#include "kernel/Coefficients/FunctionCoefficient.hpp"

#pragma once

/**
 *
 * @brief C++ function of the expression
 *
/**
 *
 * @brief C++ function of the analytical expression
 *
 * Default expression:
 *     F = (1.0/4.0)*std::pow(x1, 4) - 1.0/4.0*std::pow(x1, 2) + (1.0/4.0)*std::pow(x2, 4)
- 1.0/4.0*std::pow(x2, 2) + (1.0/4.0)*std::pow(x3, 4) - 1.0/4.0*std::pow(x3, 2)
 *
 * Conditional definitions:
 *   if (0.0 <= x1 <= 0.75)
 *       F = (1.0/4.0)*std::pow(x1, 4) - 1.0/2.0*std::pow(x1, 2) + (1.0/4.0)*std::pow(x2, 4)
- 1.0/2.0*std::pow(x2, 2) + (1.0/4.0)*std::pow(x3, 4) - 1.0/2.0*std::pow(x3, 2)
 *
 *   if (0.75 < x1 < 1)
 *       F = (1.0/4.0)*std::pow(x1, 4) - 1.0/3.0*std::pow(x1, 2) + (1.0/4.0)*std::pow(x2, 4)
- 1.0/3.0*std::pow(x2, 2) + (1.0/4.0)*std::pow(x3, 4) - 1.0/3.0*std::pow(x3, 2)
 *
 */
class ConditionalCoefficient : public FunctionCoefficient {
 private:
  double prefactor_;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const unsigned int dimension)>
  F() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  GradientF() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  HessianF() final;

 public:
  ConditionalCoefficient() { this->prefactor_ = 1.0; }
  explicit ConditionalCoefficient(const double prefactor) { this->prefactor_ = prefactor; }
  virtual ~ConditionalCoefficient() = default;
};

/**
 *
 * @brief C++ function of the expression
 *
 *   F = (1.0/4.0)*std::pow(x1, 4) - 1.0/4.0*std::pow(x1, 2) + (1.0/4.0)*std::pow(x2, 4)
 * - 1.0/4.0*std::pow(x2, 2) + (1.0/4.0)*std::pow(x3, 4) - 1.0/4.0*std::pow(x3, 2)
 *
 *   if (0.0 <= x1 <= 0.75)
 *       F = (1.0/4.0)*std::pow(x1, 4) - 1.0/2.0*std::pow(x1, 2) + (1.0/4.0)*std::pow(x2, 4)
 * - 1.0/2.0*std::pow(x2, 2) + (1.0/4.0)*std::pow(x3, 4) - 1.0/2.0*std::pow(x3, 2)
 *
 *   if (0.75 < x1 < 1)
 *       F = (1.0/4.0)*std::pow(x1, 4) - 1.0/3.0*std::pow(x1, 2) + (1.0/4.0)*std::pow(x2, 4)
 * - 1.0/3.0*std::pow(x2, 2) + (1.0/4.0)*std::pow(x3, 4) - 1.0/3.0*std::pow(x3, 2)
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const unsigned int dimension)>
ConditionalCoefficient::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x1 = input_vector[0];
    double x2 = input_vector[1];
    double x3 = input_vector[2];
    double F = (1.0 / 4.0) * std::pow(x1, 4) - 1.0 / 4.0 * std::pow(x1, 2) +
               (1.0 / 4.0) * std::pow(x2, 4) - 1.0 / 4.0 * std::pow(x2, 2) +
               (1.0 / 4.0) * std::pow(x3, 4) - 1.0 / 4.0 * std::pow(x3, 2);
    if ((0.0 <= x1) && (x1 <= 0.75)) {
      F = (1.0 / 4.0) * std::pow(x1, 4) - 1.0 / 2.0 * std::pow(x1, 2) +
          (1.0 / 4.0) * std::pow(x2, 4) - 1.0 / 2.0 * std::pow(x2, 2) +
          (1.0 / 4.0) * std::pow(x3, 4) - 1.0 / 2.0 * std::pow(x3, 2);
    }
    if ((0.75 < x1) && (x1 < 1)) {
      F = (1.0 / 4.0) * std::pow(x1, 4) - 1.0 / 3.0 * std::pow(x1, 2) +
          (1.0 / 4.0) * std::pow(x2, 4) - 1.0 / 3.0 * std::pow(x2, 2) +
          (1.0 / 4.0) * std::pow(x3, 4) - 1.0 / 3.0 * std::pow(x3, 2);
    }
    return this->prefactor_ * F;
  };
  return func;
}

/**
 *
 * @brief Gradient
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const unsigned int dimension)>
ConditionalCoefficient::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x1 = input_vector[0];
    double x2 = input_vector[1];
    double x3 = input_vector[2];
    std::vector<double> gradient(3);
    gradient[0] = this->prefactor_ * (std::pow(x1, 3) - 1.0 / 2.0 * x1);
    gradient[1] = this->prefactor_ * (std::pow(x2, 3) - 1.0 / 2.0 * x2);
    gradient[2] = this->prefactor_ * (std::pow(x3, 3) - 1.0 / 2.0 * x3);
    if ((0.0 <= x1) && (x1 <= 0.75)) {
      gradient[0] = this->prefactor_ * (std::pow(x1, 3) - x1);
      gradient[1] = this->prefactor_ * (std::pow(x2, 3) - x2);
      gradient[2] = this->prefactor_ * (std::pow(x3, 3) - x3);
    }
    if ((0.75 < x1) && (x1 < 1)) {
      gradient[0] = this->prefactor_ * (std::pow(x1, 3) - 2.0 / 3.0 * x1);
      gradient[1] = this->prefactor_ * (std::pow(x2, 3) - 2.0 / 3.0 * x2);
      gradient[2] = this->prefactor_ * (std::pow(x3, 3) - 2.0 / 3.0 * x3);
    }
    return gradient;
  };
  return func;
}

/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 *
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const unsigned int dimension)>
ConditionalCoefficient::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x1 = input_vector[0];
    double x2 = input_vector[1];
    double x3 = input_vector[2];
    std::vector<double> hessian(9);
    hessian[0] = this->prefactor_ * (3 * std::pow(x1, 2) - 1.0 / 2.0);
    hessian[1] = this->prefactor_ * (0);
    hessian[2] = this->prefactor_ * (0);
    hessian[3] = this->prefactor_ * (0);
    hessian[4] = this->prefactor_ * (3 * std::pow(x2, 2) - 1.0 / 2.0);
    hessian[5] = this->prefactor_ * (0);
    hessian[6] = this->prefactor_ * (0);
    hessian[7] = this->prefactor_ * (0);
    hessian[8] = this->prefactor_ * (3 * std::pow(x3, 2) - 1.0 / 2.0);
    if ((0.0 <= x1) && (x1 <= 0.75)) {
      hessian[0] = this->prefactor_ * (3 * std::pow(x1, 2) - 1);
      hessian[1] = this->prefactor_ * (0);
      hessian[2] = this->prefactor_ * (0);
      hessian[3] = this->prefactor_ * (0);
      hessian[4] = this->prefactor_ * (3 * std::pow(x2, 2) - 1);
      hessian[5] = this->prefactor_ * (0);
      hessian[6] = this->prefactor_ * (0);
      hessian[7] = this->prefactor_ * (0);
      hessian[8] = this->prefactor_ * (3 * std::pow(x3, 2) - 1);
    }
    if ((0.75 < x1) && (x1 < 1)) {
      hessian[0] = this->prefactor_ * (3 * std::pow(x1, 2) - 2.0 / 3.0);
      hessian[1] = this->prefactor_ * (0);
      hessian[2] = this->prefactor_ * (0);
      hessian[3] = this->prefactor_ * (0);
      hessian[4] = this->prefactor_ * (3 * std::pow(x2, 2) - 2.0 / 3.0);
      hessian[5] = this->prefactor_ * (0);
      hessian[6] = this->prefactor_ * (0);
      hessian[7] = this->prefactor_ * (0);
      hessian[8] = this->prefactor_ * (3 * std::pow(x3, 2) - 2.0 / 3.0);
    }
    return hessian;
  };
  return func;
}
