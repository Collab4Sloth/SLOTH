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

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <vector>

#include "kernel/Coefficients/FunctionCoefficient.hpp"

#pragma once

/**
 *
 * @brief Coefficient based on expression: 5.0*(x - 0.3) * (x - 0.3) * (0.7 - x) * (0.7 - x)
 *
 */
class DoubleWell : public FunctionCoefficient {
 private:
  double prefactor_;

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
  DoubleWell() { this->prefactor_ = 1.0; }
  explicit DoubleWell(const double prefactor) { this->prefactor_ = prefactor; }
  virtual ~DoubleWell() = default;
};

/**
 *
 * @brief C++ function of the expression: 5.0*(x - 0.3) * (x - 0.3) * (0.7 - x) * (0.7 - x)
 *
 * @return std::function<double(const std::vector<double>&,const std::vector<double>&)>
 */
std::function<double(const std::vector<double>&, const std::vector<double>&,
                     const unsigned int dimension)>
DoubleWell::F() {
  auto func = [&](const std::vector<double>& input_vector,
                  [[maybe_unused]] const std::vector<double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x = input_vector[0];
    double F = std::pow(7.0 / 10.0 - x, 2) * (x - 3.0 / 10.0) * (5 * x - 3.0 / 2.0);
    return this->prefactor_ * F;
  };
  return func;
}

/**
 *
 * @brief Gradient
 *
 * @return std::function<std::vector<double>(const std::vector<double>&,const std::vector<double>&,
 * const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&,
                                  const unsigned int dimension)>
DoubleWell::GradientF() {
  auto func = [&](const std::vector<double>& input_vector,
                  [[maybe_unused]] const std::vector<double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x = input_vector[0];
    std::vector<double> gradient(1);
    gradient[0] = this->prefactor_ * (5 * std::pow(7.0 / 10.0 - x, 2) * (x - 3.0 / 10.0) +
                                      std::pow(7.0 / 10.0 - x, 2) * (5 * x - 3.0 / 2.0) +
                                      (x - 3.0 / 10.0) * (2 * x - 7.0 / 5.0) * (5 * x - 3.0 / 2.0));
    return gradient;
  };
  return func;
}

/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 *
 * @return std::function<std::vector<double>(const std::vector<double>&,const std::vector<double>&,
 * const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&,
                                  const unsigned int dimension)>
DoubleWell::HessianF() {
  auto func = [&](const std::vector<double>& input_vector,
                  [[maybe_unused]] const std::vector<double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x = input_vector[0];
    std::vector<double> hessian(1);
    hessian[0] = this->prefactor_ *
                 (10 * std::pow(7.0 / 10.0 - x, 2) + 10 * (x - 3.0 / 10.0) * (2 * x - 7.0 / 5.0) +
                  2 * (x - 3.0 / 10.0) * (5 * x - 3.0 / 2.0) +
                  2 * (2 * x - 7.0 / 5.0) * (5 * x - 3.0 / 2.0));
    return hessian;
  };
  return func;
}

/**
 *
 * @brief Coefficient based on expression: dot(x,x)
 *
 */
class GradEnergy : public FunctionCoefficient {
 private:
  double prefactor_;

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
  GradEnergy() { this->prefactor_ = 1.0; }
  explicit GradEnergy(const double prefactor) { this->prefactor_ = prefactor; }
  virtual ~GradEnergy() = default;
};

/**
 *
 * @brief C++ function of the expression: dot(x,x)
 *
 * @return std::function<double(const std::vector<double>&,const std::vector<double>&)>
 */
std::function<double(const std::vector<double>&, const std::vector<double>&,
                     const unsigned int dimension)>
GradEnergy::F() {
  auto func = [&](const std::vector<double>& input_vector,
                  [[maybe_unused]] const std::vector<double>&, const unsigned int dimension) {
    std::vector<double> x;
    for (unsigned int i = 0; i < dimension; i++) x.push_back(input_vector[0 * dimension + i]);
    double F = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    return this->prefactor_ * F;
  };
  return func;
}

/**
 *
 * @brief Gradient
 *
 * @return std::function<std::vector<double>(const std::vector<double>&,const std::vector<double>&,
 * const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&,
                                  const unsigned int dimension)>
GradEnergy::GradientF() {
  auto func = [&]([[maybe_unused]] const std::vector<double>& input_vector,
                  [[maybe_unused]] const std::vector<double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> gradient(1, 0.0);
    return gradient;
  };
  return func;
}

/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 *
 * @return std::function<std::vector<double>(const std::vector<double>&,const std::vector<double>&,
 * const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&,
                                  const unsigned int dimension)>
GradEnergy::HessianF() {
  auto func = [&]([[maybe_unused]] const std::vector<double>& input_vector,
                  [[maybe_unused]] const std::vector<double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> hessian(1, 0.0);
    return hessian;
  };
  return func;
}
