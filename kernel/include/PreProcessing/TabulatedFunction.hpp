/**
 * @file TabulatedFunction.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Reads a CSV table (t, x1..xN, F) sampled on a structured grid
 *        and provides fast (multi)linear interpolation on it.
 *        The number of inputs N (1, 2 or 3) is auto-detected from the CSV.
 *        No geometry conversion is applied: points are given in CSV coordinates.
 * @version 0.1
 * @date 2026-09-22
 *
 * @copyright CEA (C) 2026
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
#include <array>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "Options/Options.hpp"
#include "mfem.hpp"  // NOLINT

/**
 * @brief Reads a CSV table (t, x1..xN, F) sampled on a structured grid
 *        and provides fast (multi)linear interpolation on it.
 *        The number of inputs N (1, 2 or 3) is auto-detected from the CSV.
 *        No geometry conversion is applied: points are given in CSV coordinates.
 */
class TabulatedFunction {
 private:
  static constexpr int kMaxDim = 3;

  struct Sample {
    double t;
    std::array<double, kMaxDim> x{};
    double F;
  };

  std::vector<Sample> read_csv(const std::filesystem::path& path);
  void build_grid(std::vector<Sample>& samples);

  static std::vector<double> unique_sorted(std::vector<double>& values);
  std::pair<int, double> locate(const std::vector<double>& axis, bool uniform, double val) const;
  std::size_t flat_index(int it, const std::array<int, kMaxDim>& idx) const;

  int spatial_dim_ = 0;
  std::vector<double> t_axis_;
  bool t_uniform_ = false;
  std::array<std::vector<double>, kMaxDim> x_axis_;
  std::array<bool, kMaxDim> x_uniform_{};
  std::vector<double> data_;

 public:
  explicit TabulatedFunction(const std::string& csv_path);
  double interpolate(double t, const mfem::Vector& pos) const;
  double interpolate(double t, std::span<const double> pos) const;

  int input_dim() const { return spatial_dim_; }
  bool has_data() const { return !data_.empty(); }
};
