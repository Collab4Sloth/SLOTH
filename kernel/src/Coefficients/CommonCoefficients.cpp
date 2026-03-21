/**
 * @file CommonCoefficients.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief
 * @version 0.1
 * @date 2026-03-21
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

#include "Coefficients/CommonCoefficients.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <span>
#include <vector>

#include "Coefficients/FunctionCoefficient.hpp"

/**
 *
 * @brief C++ function of the expression: x*x*(1.0-x) * (1.0-x)
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const unsigned int dimension)>
W::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x = input_vector[0];
    double F = std::pow(x, 2) * std::pow(1 - x, 2);
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
W::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x = input_vector[0];
    std::vector<double> gradient(1);
    gradient[0] = this->prefactor_ * (std::pow(x, 2) * (2 * x - 2) + 2 * x * std::pow(1 - x, 2));
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
W::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x = input_vector[0];
    std::vector<double> hessian(1);
    hessian[0] =
        this->prefactor_ * (2 * std::pow(x, 2) + 4 * x * (2 * x - 2) + 2 * std::pow(1 - x, 2));
    return hessian;
  };
  return func;
}

/**
 *
 * @brief C++ function of the expression: 0.25 * (x * x - 1.0) * (x * x - 1.0)
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const unsigned int dimension)>
Fw::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x = input_vector[0];
    double F = ((1.0 / 4.0) * std::pow(x, 2) - 1.0 / 4.0) * (std::pow(x, 2) - 1);
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
Fw::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x = input_vector[0];
    std::vector<double> gradient(1);
    gradient[0] = this->prefactor_ * (2 * x * ((1.0 / 4.0) * std::pow(x, 2) - 1.0 / 4.0) +
                                      (1.0 / 2.0) * x * (std::pow(x, 2) - 1));
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
Fw::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x = input_vector[0];
    std::vector<double> hessian(1);
    hessian[0] = this->prefactor_ * (3 * std::pow(x, 2) - 1);
    return hessian;
  };
  return func;
}

/**
 *
 * @brief C++ function of the expression: x * x * x * (6.0 * x * x - 15.0 * x + 10.0)
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const unsigned int dimension)>
H::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x = input_vector[0];
    double F = std::pow(x, 3) * (6 * std::pow(x, 2) - 15 * x + 10);
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
H::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x = input_vector[0];
    std::vector<double> gradient(1);
    gradient[0] = this->prefactor_ * (std::pow(x, 3) * (12 * x - 15) +
                                      3 * std::pow(x, 2) * (6 * std::pow(x, 2) - 15 * x + 10));
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
H::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double x = input_vector[0];
    std::vector<double> hessian(1);
    hessian[0] = this->prefactor_ * (12 * std::pow(x, 3) + 6 * std::pow(x, 2) * (12 * x - 15) +
                                     6 * x * (6 * std::pow(x, 2) - 15 * x + 10));
    return hessian;
  };
  return func;
}

/**
 *
 * @brief C++ function of the expression: x * log(x) + (1.0 - x)*log(1.0 - x)
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const unsigned int dimension)>
Log::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double epsilon = 1.e-10;
    double x = std::min(1.0 - epsilon, std::max(epsilon, input_vector[0]));
    double F = x * std::log(x) + (1 - x) * std::log(1 - x);
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
Log::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double epsilon = 1.e-10;
    double x = std::min(1.0 - epsilon, std::max(epsilon, input_vector[0]));
    std::vector<double> gradient(1);
    gradient[0] = this->prefactor_ * (std::log(x) - std::log(1 - x));
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
Log::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    double epsilon = 1.e-10;
    double x = std::min(1.0 - epsilon, std::max(epsilon, input_vector[0]));
    std::vector<double> hessian(1);
    hessian[0] = this->prefactor_ * (1.0 / (1 - x) + 1.0 / x);
    return hessian;
  };
  return func;
}

/**
 *
 * @brief C++ function of the expression: 0.5*dot(x,x)
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&)>
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const unsigned int dimension)>
GradientEnergy::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&, const unsigned int dimension) {
    std::vector<double> x;
    for (unsigned int i = 0; i < dimension; i++) x.push_back(input_vector[0 * dimension + i]);
    double F = 0.5 * std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
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
GradientEnergy::GradientF() {
  auto func = [&]([[maybe_unused]] const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
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
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const
 * double>&, const unsigned int dimension)>
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const unsigned int dimension)>
GradientEnergy::HessianF() {
  auto func = [&]([[maybe_unused]] const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> hessian(1, 0.0);
    return hessian;
  };
  return func;
}
