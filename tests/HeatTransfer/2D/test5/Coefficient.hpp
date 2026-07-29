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
 *       F = -pi*t*cos(pi*x)
 */
class NeumannCoefficient : public FunctionCoefficient {
 private:
  double prefactor_;
 protected:
  std::function<double(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> F() final;
  std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> GradientF() final;
  std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&,const std::span<const double>&, const unsigned int dimension)> HessianF() final;

 public:
  NeumannCoefficient() : prefactor_(1.0) {}
 explicit NeumannCoefficient(const double prefactor): prefactor_(prefactor) {}
  virtual ~NeumannCoefficient() = default;
};

/**
 *
 * @brief C++ function of the expression
 *
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)>
 */
 std::function<double(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> NeumannCoefficient::F() {
  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&,const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {
    double T = input_vector[0];
    double x = auxiliary_vector[0];
    double y = auxiliary_vector[1];
    double t = this->time_;
 double   F = -M_PI*t*std::cos(M_PI*x);
    return this->prefactor_ * F;
  };
  return func;
}

/**
 *
 * @brief Gradient
 * 
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const unsigned int dimension)> 
 */
std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> NeumannCoefficient::GradientF() {
  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&,const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {
    double T = input_vector[0];
    double x = auxiliary_vector[0];
    double y = auxiliary_vector[1];
    double t = this->time_;
    std::vector<double> gradient(1);
    gradient[0] = this->prefactor_ * (0);
    return gradient;
  };
  return func;
}

/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 * 
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> 
 */
std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> NeumannCoefficient::HessianF() {
  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&,const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {
    double T = input_vector[0];
    double x = auxiliary_vector[0];
    double y = auxiliary_vector[1];
    double t = this->time_;
    std::vector<double> hessian(1);
    hessian[0] = this->prefactor_ * (0);
    return hessian;
  };
  return func;
}
/**
 *
 * @brief C++ function of the analytical expression
 *
 *       F = t*sin(pi*y)
 */
class DirichletCoefficient : public FunctionCoefficient {
 private:
  double prefactor_;
 protected:
  std::function<double(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> F() final;
  std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> GradientF() final;
  std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&,const std::span<const double>&, const unsigned int dimension)> HessianF() final;

 public:
  DirichletCoefficient() : prefactor_(1.0) {}
 explicit DirichletCoefficient(const double prefactor): prefactor_(prefactor) {}
  virtual ~DirichletCoefficient() = default;
};

/**
 *
 * @brief C++ function of the expression
 *
 *
 * @return std::function<double(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)>
 */
 std::function<double(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> DirichletCoefficient::F() {
  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&,const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {
    double T = input_vector[0];
    double x = auxiliary_vector[0];
    double y = auxiliary_vector[1];
    double t = this->time_;
 double   F = t*std::sin(M_PI*y);
    return this->prefactor_ * F;
  };
  return func;
}

/**
 *
 * @brief Gradient
 * 
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const unsigned int dimension)> 
 */
std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> DirichletCoefficient::GradientF() {
  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&,const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {
    double T = input_vector[0];
    double x = auxiliary_vector[0];
    double y = auxiliary_vector[1];
    double t = this->time_;
    std::vector<double> gradient(1);
    gradient[0] = this->prefactor_ * (0);
    return gradient;
  };
  return func;
}

/**
 *
 * @brief Hessian
 * @remark Hessian matrix stored in vector : H(i,j)->H(i*n+j)
 * 
 * @return std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> 
 */
std::function<std::vector<double>(const std::span<const double>&,const std::span<const double>&, const std::span<const double>&, const unsigned int dimension)> DirichletCoefficient::HessianF() {
  auto func = [&](const std::span<const double>& input_vector, [[maybe_unused]] const std::span<const double>&,const std::span<const double>& auxiliary_vector, [[maybe_unused]] const unsigned int dimension) {
    double T = input_vector[0];
    double x = auxiliary_vector[0];
    double y = auxiliary_vector[1];
    double t = this->time_;
    std::vector<double> hessian(1);
    hessian[0] = this->prefactor_ * (0);
    return hessian;
  };
  return func;
}
