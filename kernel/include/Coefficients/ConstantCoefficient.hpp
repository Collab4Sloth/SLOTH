/**
 * @file ConstantCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Constant coefficient
 * @version 0.1
 * @date 2025-09-05
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
#include <cmath>
#include <functional>
#include <span>
#include <vector>

#include "Coefficients/FunctionCoefficient.hpp"

/**
 * @class ConstantCoefficient
 * @brief Constant-valued coefficient implementation.
 *
 * This class represents a constant scalar-valued coefficient.
 *
 * The class derives from FunctionCoefficient and overrides the
 * function, gradient, and Hessian evaluators.
 */
class ConstantCoefficient : public FunctionCoefficient {
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
  double value_;
  inline explicit ConstantCoefficient(const double value) : value_(value) {}
  ~ConstantCoefficient() override = default;
};
