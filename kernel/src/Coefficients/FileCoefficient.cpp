/**
 * @file FileCoefficient.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class defining a FunctionCoefficient from file
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

#include "Coefficients/FileCoefficient.hpp"

#include <any>
#include <array>
#include <cmath>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/FunctionCoefficient.hpp"
#include "Options/Options.hpp"
#include "PreProcessing/TabulatedFunction.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new FileCoefficient::FileCoefficient object
 *
 * @param filename CSV file
 * @param geom Geometry (Cartesian or Axisymmetric)
 * @param axis_center axes of symmetry
 */
FileCoefficient::FileCoefficient(const std::string& filename, Geometry geom,
                                 std::array<double, 2> axis_center)
    : prefactor_(1.0), geometry_(geom), center_(axis_center) {
  this->file_coefficient_ = std::filesystem::path(filename);
  std::string msg = "Error while trying to open " + filename;
  MFEM_VERIFY(std::filesystem::exists(this->file_coefficient_), msg.c_str());
  this->field_ = std::make_shared<TabulatedFunction>(this->file_coefficient_);
  MFEM_VERIFY(this->geometry_ != Geometry::Axisymmetric || this->field_->input_dim() <= 2,
              "FileCoefficient: axisymmetric CSV must be (t, r, F) or (t, r, z, F).");
}
/**
 * @brief Geometry of the tabulated field.
 *
 * @return Geometry Cartesian, or Axisymmetric if mesh points are mapped to
 *         r = sqrt((x - xc)^2 + (y - yc)^2) before interpolation.
 */
Geometry FileCoefficient::geometry() const { return geometry_; }

/**
 * @brief Position of the symmetry axis in the mesh frame.
 *
 * @return const std::array<double, 2>& Coordinates (xc, yc) of the axis,
 *         normal to the (x, y) plane. Meaningful only in axisymmetric mode.
 */
const std::array<double, 2>& FileCoefficient::axis_center() const { return center_; }

/**
 * @brief Number of input columns in the CSV file.
 *
 * Number of coordinates between the time and the value columns, i.e.
 * N for a file `t, x1, ..., xN, F` (1 for `t, r, F`, 2 for `t, r, z, F`).
 *
 * @return int Number of input coordinates of the tabulated function.
 */
int FileCoefficient::input_dim() const { return field_->input_dim(); }

/**
 * @brief Value of the coefficient as a FunctionCoefficient function.
 *
 * The first argument (input_vector) must hold the physical mesh point; the
 * evaluation time is the one set with set_time().
 */
std::function<double(const std::span<const double>&, const std::span<const double>&,
                     const std::span<const double>&, const unsigned int dimension)>
FileCoefficient::F() {
  auto func = [this](const std::span<const double>& input_vector,
                     [[maybe_unused]] const std::span<const double>&,
                     [[maybe_unused]] const std::span<const double>&,
                     [[maybe_unused]] const unsigned int dimension) {
    return this->value(this->time_, input_vector);
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
FileCoefficient::GradientF() {
  return [](const std::span<const double>&, const std::span<const double>&,
            const std::span<const double>&, const unsigned int) -> std::vector<double> {
    MFEM_ABORT("FileCoefficient: gradient is not available for a coefficient read from file.");
    return {};
  };
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
FileCoefficient::HessianF() {
  return [](const std::span<const double>&, const std::span<const double>&,
            const std::span<const double>&, const unsigned int) -> std::vector<double> {
    MFEM_ABORT("FileCoefficient: Hessian is not available for a coefficient read from file.");
    return {};
  };
}

/**
 * @brief Set time for all FunctionCoefficients
 *
 * @param time
 */
void FileCoefficient::set_time(double time) { this->time_ = time; }

/**
 * @brief Evaluate the coefficient at a given time and physical mesh point.
 *
 * Maps the mesh point to the CSV input coordinates according to the geometry,
 * then interpolates the tabulated function:
 *  - Cartesian: coordinates are used as is;
 *  - Axisymmetric: r = sqrt((x - xc)^2 + (y - yc)^2), and z = coordinates[2]
 *    is appended when the CSV is (t, r, z, F).
 *
 * @param time        Evaluation time.
 * @param coordinates Physical point in the mesh frame.
 * @return double prefactor * F(time, mapped coordinates).
 */
double FileCoefficient::value(double time, std::span<const double> coordinates) const {
  if (this->geometry_ == Geometry::Cartesian) {
    return this->prefactor_ * this->field_->interpolate(time, coordinates);
  }

  const int n = this->field_->input_dim();
  MFEM_VERIFY(static_cast<int>(coordinates.size()) >= n + 1,
              "FileCoefficient: mesh dimension too small for this axisymmetric field ("
                  << coordinates.size() << " coordinates, CSV input dimension = " << n << ")");
  std::array<double, 2> q{
      std::hypot(coordinates[0] - this->center_[0], coordinates[1] - this->center_[1]), 0.0};
  if (n == 2) q[1] = coordinates[2];
  return this->prefactor_ * this->field_->interpolate(time, std::span<const double>(q.data(), n));
}