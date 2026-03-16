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

#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <span>
#include <vector>

#include "Coefficients/FunctionCoefficient.hpp"

/**
 *
 * @brief Coefficient based on expression: x*x*(1.0-x) * (1.0-x)
 *
 */
class W : public FunctionCoefficient {
 private:
  double prefactor_;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const unsigned int dimension)>
  F() override final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  GradientF() override final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  HessianF() override final;

 public:
  W() { this->prefactor_ = 1.0; }
  explicit W(const double prefactor) { this->prefactor_ = prefactor; }
  virtual ~W() = default;
};

/**
 *
 * @brief Coefficient based on expression: 0.25 * (x * x - 1.0) * (x * x - 1.0)
 *
 */
class Fw : public FunctionCoefficient {
 private:
  double prefactor_;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const unsigned int dimension)>
  F() override final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  GradientF() override final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  HessianF() override final;

 public:
  Fw() { this->prefactor_ = 1.0; }
  explicit Fw(const double prefactor) { this->prefactor_ = prefactor; }
  virtual ~Fw() = default;
};

/**
 *
 * @brief Coefficient based on expression: x * x * x * (6.0 * x * x - 15.0 * x + 10.0)
 *
 */
class H : public FunctionCoefficient {
 private:
  double prefactor_;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const unsigned int dimension)>
  F() override final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  GradientF() override final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  HessianF() override final;

 public:
  H() { this->prefactor_ = 1.0; }
  explicit H(const double prefactor) { this->prefactor_ = prefactor; }
  virtual ~H() = default;
};

/**
 *
 * @brief Coefficient based on expression: x * log(x) + (1.0 - x)*log(1.0 - x)
 *
 */
class Log : public FunctionCoefficient {
 private:
  double prefactor_;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const unsigned int dimension)>
  F() override final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  GradientF() override final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  HessianF() override final;

 public:
  Log() { this->prefactor_ = 1.0; }
  explicit Log(const double prefactor) { this->prefactor_ = prefactor; }
  virtual ~Log() = default;
};

/**
 *
 * @brief Coefficient based on expression: 0.5*dot(x,x)
 *
 */
class GradientEnergy : public FunctionCoefficient {
 private:
  double prefactor_;

 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const unsigned int dimension)>
  F() override final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  GradientF() override final;
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const unsigned int dimension)>
  HessianF() override final;

 public:
  GradientEnergy() { this->prefactor_ = 1.0; }
  explicit GradientEnergy(const double prefactor) { this->prefactor_ = prefactor; }
  virtual ~GradientEnergy() = default;
};
