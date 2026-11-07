/**
 * @file CommonCoefficients.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
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
 * @brief Value of W (see WEval)
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
W::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double value = 0.0;
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, &value, nullptr, nullptr);
    return value;
  };
  return func;
}

/**
 * @brief Gradient of W (see WEval)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
W::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> gradient(1);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, gradient.data(), nullptr);
    return gradient;
  };
  return func;
}

/**
 * @brief Hessian of W (see WEval)
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
W::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> hessian(1);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, nullptr, hessian.data());
    return hessian;
  };
  return func;
}

/**
 * @brief Value of Fw (see FwEval)
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
Fw::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double value = 0.0;
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, &value, nullptr, nullptr);
    return value;
  };
  return func;
}

/**
 * @brief Gradient of Fw (see FwEval)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Fw::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> gradient(1);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, gradient.data(), nullptr);
    return gradient;
  };
  return func;
}

/**
 * @brief Hessian of Fw (see FwEval)
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Fw::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> hessian(1);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, nullptr, hessian.data());
    return hessian;
  };
  return func;
}

/**
 * @brief Value of H (see HEval)
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
H::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double value = 0.0;
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, &value, nullptr, nullptr);
    return value;
  };
  return func;
}

/**
 * @brief Gradient of H (see HEval)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
H::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> gradient(1);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, gradient.data(), nullptr);
    return gradient;
  };
  return func;
}

/**
 * @brief Hessian of H (see HEval)
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
H::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> hessian(1);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, nullptr, hessian.data());
    return hessian;
  };
  return func;
}

/**
 * @brief Value of Log (see LogEval)
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
Log::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    double value = 0.0;
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, &value, nullptr, nullptr);
    return value;
  };
  return func;
}

/**
 * @brief Gradient of Log (see LogEval)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Log::GradientF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> gradient(1);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, gradient.data(), nullptr);
    return gradient;
  };
  return func;
}

/**
 * @brief Hessian of Log (see LogEval)
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 */
std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                  const std::span<const double>&, const unsigned int dimension)>
Log::HessianF() {
  auto func = [&](const std::span<const double>& input_vector,
                  const std::span<const double>& exp_input_vector,
                  const std::span<const double>& auxiliary_vector,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> hessian(1);
    this->eval_(input_vector.data(), exp_input_vector.data(), auxiliary_vector.data(),
                this->coordinates_.data(), this->time_, nullptr, nullptr, hessian.data());
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
                     const std::span<const double>&, const unsigned int dimension)>
GradientEnergy::F() {
  auto func = [&](const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
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
                                  const std::span<const double>&, const unsigned int dimension)>
GradientEnergy::GradientF() {
  auto func = [&]([[maybe_unused]] const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
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
                                  const std::span<const double>&, const unsigned int dimension)>
GradientEnergy::HessianF() {
  auto func = [&]([[maybe_unused]] const std::span<const double>& input_vector,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const std::span<const double>&,
                  [[maybe_unused]] const unsigned int dimension) {
    std::vector<double> hessian(1, 0.0);
    return hessian;
  };
  return func;
}
