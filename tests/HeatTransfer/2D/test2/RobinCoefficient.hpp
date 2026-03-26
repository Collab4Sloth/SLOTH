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

#include "Coefficients/FunctionCoefficient.hpp"

#pragma once

/**
 *
 * @brief C++ function of the analytical expression
 *
 *       F = T**3*epsilon*sigma + h
 */
class RobinCoefficient : public FunctionCoefficient {
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
  RobinCoefficient() : prefactor_(1.0) {}
  explicit RobinCoefficient(const double prefactor) : prefactor_(prefactor) {}
  virtual ~RobinCoefficient() = default;
};

/**
 *
 * @brief C++ function of the expression
 *
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const unsigned int dimension)>
RobinCoefficient::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double T = input_vector[0];
    double h = 50.0;
    double epsilon = 0.7;
    double sigma = 5.669e-8;
    double F = std::pow(T, 3) * epsilon * sigma + h;
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
RobinCoefficient::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double T = input_vector[0];
    double h = 50.0;
    double epsilon = 0.7;
    double sigma = 5.669e-8;
    std::vector<double> gradient(1);
    gradient[0] = this->prefactor_ * (3 * std::pow(T, 2) * epsilon * sigma);
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
RobinCoefficient::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double T = input_vector[0];
    double h = 50.0;
    double epsilon = 0.7;
    double sigma = 5.669e-8;
    std::vector<double> hessian(1);
    hessian[0] = this->prefactor_ * (6 * T * epsilon * sigma);
    return hessian;
  };
  return func;
}
