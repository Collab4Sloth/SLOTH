/**
 * @file ConstantCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Constant coefficient
 * @version 0.1
 * @date 2025-09-05
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
 * @brief Constant Coefficient
 *
 */
class ConstantCoefficient : public FunctionCoefficient {
 protected:
  std::function<double(const std::vector<double>&, const std::vector<double>&,
                       const unsigned int dimension)>
  F() final;
  std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&,
                                    const unsigned int dimension)>
  GradientF() final;
  std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&,
                                    const unsigned int dimension)>
  HessianF() final;

 public:
  double value_;
  inline explicit ConstantCoefficient(const double value) : value_(value) {}
  ~ConstantCoefficient() override = default;
};

/**
 *
 * @brief C++ function of the expression: constant
 *
 * @return std::function<double(const std::vector<double>&)>
 */
std::function<double(const std::vector<double>&, const std::vector<double>&,
                     const unsigned int dimension)>
ConstantCoefficient::F() {
  auto func = [&]([[maybe_unused]] const std::vector<double>& input_vector,
                  [[maybe_unused]] const std::vector<double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    // CCI
    return this->value_;
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
                                  const unsigned int dimension)>
ConstantCoefficient::GradientF() {
  auto func = [](const std::vector<double>& input_vector,
                 [[maybe_unused]] const std::vector<double>&,
                 [[maybe_unused]] const unsigned int dimension) {
    const size_t size = input_vector.size();
    std::vector<double> gradient(size, 0.0);
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
                                  const unsigned int dimension)>
ConstantCoefficient::HessianF() {
  auto func = [](const std::vector<double>& input_vector,
                 [[maybe_unused]] const std::vector<double>&,
                 [[maybe_unused]] const unsigned int dimension) {
    const size_t size = input_vector.size();
    std::vector<double> hessian(size, 0.0);
    return hessian;
  };
  return func;
}
