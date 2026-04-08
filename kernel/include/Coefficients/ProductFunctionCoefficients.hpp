/**
 * @file ProductFunctionCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class defining analytical expression, first and second derivatives of a Product of
 * FunctionCoefficient
 * @version 0.1
 * @date 2025-02-03
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
#pragma once
#include <any>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/FunctionCoefficient.hpp"

/**
 * @brief Define a product of FunctionCoefficient objects as a FunctionCoefficient
 * @remark This suggest that all FunctionCoefficient share the same variables...
 *
 */
class ProductCoefficient : public FunctionCoefficient {
 private:
  std::vector<FunctionCoefficient*> vect_coefficients_;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int dimension)>
  F() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  GradientF() final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int dimension)>
  HessianF() final;

 public:
  explicit ProductCoefficient(std::initializer_list<FunctionCoefficient*> coeffs);
  ProductCoefficient() = default;

  virtual ~ProductCoefficient() = default;

  void set_time(double time) override;
};
