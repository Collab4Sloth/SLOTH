/**
 * @file FileCoefficient.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class defining FunctionCoefficient from file
 * @version 0.1
 * @date 2026-09-22
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
#pragma once
#include <any>
#include <cmath>
#include <filesystem>
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
#include "Options/Options.hpp"
#include "PreProcessing/TabulatedFunction.hpp"

/**
 * @brief FunctionCoefficient defined by a tabulated function read from a CSV file.
 *
 * The coefficient is evaluated at a physical mesh point, mapped to the CSV
 * input coordinates according to the geometry (Cartesian or Axisymmetric).
 * Through ::Coefficient::compute(), the mesh point must be passed as the first
 * argument (input_vector), as done by MfemFileCoefficient.
 *
 * @note Only the value is available: gradient and Hessian abort if evaluated.
 */
class FileCoefficient : public FunctionCoefficient {
 private:
  double prefactor_ = 1.0;
  std::filesystem::path file_coefficient_;
  std::shared_ptr<TabulatedFunction> field_;

  Geometry geometry_ = Geometry::Cartesian;
  std::array<double, 2> center_{0.0, 0.0};

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
  FileCoefficient(const std::string& filename, Geometry geom = Geometry::Cartesian,
                  std::array<double, 2> axis_center = {0.0, 0.0});

  Geometry geometry() const;
  const std::array<double, 2>& axis_center() const;
  int input_dim() const;

  double value(double time, std::span<const double> coordinates) const;
  virtual ~FileCoefficient() = default;

  void set_time(double time) override;
};
