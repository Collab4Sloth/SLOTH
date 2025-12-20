/**
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

#include <cmath>
#include <functional>
#include <vector>

#include "kernel/Coefficients/FunctionCoefficient.hpp"

#pragma once

/**
 *
 * @brief Coefficient based on expression: x*log(x)+(1.0-x)*log(1.0-x)
 *
 */
class energylog : public FunctionCoefficient {
 protected:
  std::function<double(const std::vector<double>&, const std::vector<double>&, const int dimension)>
  F() final;
  std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&,
                                    const int dimension)>
  GradientF() final;
  std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&,
                                    const int dimension)>
  HessianF() final;

 public:
  energylog() = default;
  ~energylog() override = default;
};

/**
 *
 * @brief C++ function of the expression: x*log(x)+(1.0-x)*log(1.0-x)
 *
 * @return std::function<double(const std::vector<double>&)>
 */
std::function<double(const std::vector<double>&, const std::vector<double>&, const int dimension)>
energylog::F() {
  auto func = [](const std::vector<double>& input_vector,
                 [[maybe_unused]] const std::vector<double>&,
                 [[maybe_unused]] const int dimension) {
    double x = input_vector[0];
    double F = x * std::log(x) + (1 - x) * std::log(1 - x);
    return F;
  };
  return func;
}

/**
 *
 * @brief Gradient
 *
 * @return std::function<std::vector<double>(const std::vector<double>&)>
 */
std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&,
                                  const int dimension)>
energylog::GradientF() {
  auto func = [](const std::vector<double>& input_vector,
                 [[maybe_unused]] const std::vector<double>&,
                 [[maybe_unused]] const int dimension) {
    double x = input_vector[0];
    std::vector<double> gradient(1);
    gradient[0] = std::log(x) - std::log(1 - x);
    return gradient;
  };
  return func;
}

/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 *
 * @return std::function<std::vector<double>(const std::vector<double>&)>
 */
std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&,
                                  const int dimension)>
energylog::HessianF() {
  auto func = [](const std::vector<double>& input_vector,
                 [[maybe_unused]] const std::vector<double>&,
                 [[maybe_unused]] const int dimension) {
    double x = input_vector[0];
    std::vector<double> hessian(1);
    hessian[0] = 1.0 / (1 - x) + 1.0 / x;
    return hessian;
  };
  return func;
}
