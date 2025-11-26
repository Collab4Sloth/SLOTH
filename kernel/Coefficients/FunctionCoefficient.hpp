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
#pragma once
class FunctionCoefficient {
 protected:
  virtual std::function<double(const std::vector<double>&)> F() = 0;
  virtual std::function<std::vector<double>(const std::vector<double>&)> GradientF() = 0;
  virtual std::function<std::vector<double>(const std::vector<double>&)> HessianF() = 0;

 public:
  FunctionCoefficient();
  ~FunctionCoefficient();

  double eval_f(const std::vector<double>& input_vector);
  double eval_gradient(const int i, const std::vector<double>& input_vector);
  double eval_hessian(const int i, const int j, const std::vector<double>& input_vector);
};

FunctionCoefficient::FunctionCoefficient() {}

FunctionCoefficient::~FunctionCoefficient() {}
/**
 * @brief Compute function
 *
 * @param input_vector
 * @return double
 */
double FunctionCoefficient::eval_f(const std::vector<double>& input_vector) {
  return F()(input_vector);
}
/**
 * @brief Compute gradient at index i
 *
 * @param i
 * @param input_vector
 * @return double
 */
double FunctionCoefficient::eval_gradient(const int i, const std::vector<double>& input_vector) {
  return GradientF()(input_vector)[i];
}
/**
 * @brief Compute the Hessian at position (i,j)
 *
 * @param i
 * @param j
 * @param input_vector
 * @return double
 */
double FunctionCoefficient::eval_hessian(const int i, const int j,
                                         const std::vector<double>& input_vector) {
  const int size = input_vector.size();
  return HessianF()(input_vector)[i * size + j];
}
