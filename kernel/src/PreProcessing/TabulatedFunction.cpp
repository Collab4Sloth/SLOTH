/**
 * @file TabulatedFunction.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Reads a CSV table (t, x1..xN, F) sampled on a structured grid
 *         and provides fast (multi)linear interpolation on it.
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

#include "PreProcessing/TabulatedFunction.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <sstream>

#include "Options/Options.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a tabulated function from a CSV file.
 *
 * Reads the samples from csv_path with read_csv(), then reorganizes them into
 * a rectilinear (time × inputs) grid with build_grid(), which sets up the
 * axes, detects uniform spacing and fills data_.
 *
 * @param csv_path    Path to the CSV file containing the samples.
 *
 * @throws std::runtime_error if the CSV is malformed or does not describe a
 *         complete structured grid (see read_csv() and build_grid()).
 */
TabulatedFunction::TabulatedFunction(const std::string& csv_path) {
  std::vector<Sample> samples = read_csv(csv_path);
  build_grid(samples);
}

/**
 * @brief Read the samples of a tabulated function from a CSV file.
 *
 * Expected format:
 *  - first line: header (skipped, content not checked);
 *  - following lines: comma-separated values `t, x1, ..., xN, F`, where t is
 *    the time, x1..xN the input coordinates and F the function value;
 *  - the time column is mandatory. For a stationary field, use a single
 *    time value (e.g. t = 0 on every line): the value is then returned
 *    unchanged whatever the evaluation time;
 *  - empty lines are ignored.
 *
 * The number of inputs N is deduced from the number of columns in the first
 * data line (N = columns - 2), and all subsequent lines must have the same
 * number of columns.
 *
 * @param path Path to the CSV file.
 * @return std::vector<TabulatedFunction::Sample> Samples in file order
 *         (unsorted, not deduplicated).
 *
 * @post spatial_dim_ is set from the column count (side effect on the object).
 *
 * @throws std::runtime_error if a data line has a different number of columns
 *         than the first one.
 * @throws std::invalid_argument / std::out_of_range (from std::stod) if a value
 *         is not a valid number.
 * @note Aborts via MFEM_VERIFY if the file cannot be opened or if the deduced
 *       number of inputs is not in [1, kMaxDim].
 */
std::vector<TabulatedFunction::Sample> TabulatedFunction::read_csv(
    const std::filesystem::path& path) {
  std::ifstream file(path);
  MFEM_VERIFY(file.is_open(), "TabulatedFunction: cannot open CSV file '" << path.string() << "'.");

  std::vector<Sample> samples;
  std::string line;
  // Read header
  std::getline(file, line);

  int expectedColumns = -1;
  int lineNo = 1;

  while (std::getline(file, line)) {
    ++lineNo;
    if (line.empty()) continue;

    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) fields.push_back(token);

    if (expectedColumns < 0) {
      expectedColumns = static_cast<int>(fields.size());
      spatial_dim_ = expectedColumns - 2;
      MFEM_VERIFY(spatial_dim_ >= 1 && spatial_dim_ <= kMaxDim,
                  "TabulatedFunction: input coordinates not well read in CSV file. Please check "
                  "data.");
    } else if (static_cast<int>(fields.size()) != expectedColumns) {
      throw std::runtime_error("TabulatedFunction: line " + std::to_string(lineNo) + " contains " +
                               std::to_string(fields.size()) + " columns, " +
                               std::to_string(expectedColumns) + " expected");
    }

    Sample s;
    s.t = std::stod(fields[0]);
    for (int d = 0; d < spatial_dim_; ++d) s.x[d] = std::stod(fields[1 + d]);
    s.F = std::stod(fields[1 + spatial_dim_]);
    samples.push_back(s);
  }
  MFEM_VERIFY(!samples.empty(), "TabulatedFunction: CSV file '"
                                    << path.string()
                                    << "' contains no data line (only a header, or empty file).");
  return samples;
}

/**
 * @brief Sort values and remove exact duplicates.
 *
 * Sorts values in ascending order and removes consecutive duplicates with
 * std::unique, yielding the strictly increasing set of distinct values.
 * Used to extract the axis coordinates from the samples.
 *
 * @param values Values to process. Modified in place: on return it holds the
 *               sorted, deduplicated values.
 * @return std::vector<double> A copy of the sorted, deduplicated values.
 *
 * @note Duplicates are detected by exact floating-point equality. Coordinates
 *       that differ only by round-off (e.g. 0.1 vs 0.1000000000000001) are
 *       kept as distinct values.
 */
std::vector<double> TabulatedFunction::unique_sorted(std::vector<double>& values) {
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
  return values;
}

/**
 * @brief Build the rectilinear (time × inputs) grid from the samples.
 *
 * Extracts the distinct time values and the distinct coordinates along each
 * input dimension to form t_axis_ and x_axis_[d] (sorted, strictly
 * increasing), detects whether each axis is uniformly spaced (t_uniform_,
 * x_uniform_[d]), then fills data_ in the layout expected by flat_index()
 * (time slowest, last input dimension fastest).
 *
 * The grid is accepted only if every node receives exactly one sample:
 * the number of samples must equal the number of nodes, and no
 * (t, x1, ..., xN) tuple may appear twice.
 *
 * @param samples Samples read from the CSV file.
 *
 * @pre samples is non-empty and all samples share the same number of inputs.
 * @post t_axis_, x_axis_, t_uniform_, x_uniform_ and data_ are initialized
 *       and consistent with flat_index().
 *
 * @throws std::runtime_error if the number of samples does not match the
 *         number of grid nodes (missing points or non-structured data).
 * @throws std::runtime_error if a (t, x1, ..., xN) tuple appears twice.
 */
void TabulatedFunction::build_grid(std::vector<Sample>& samples) {
  std::vector<double> tRaw;
  tRaw.reserve(samples.size());
  std::array<std::vector<double>, kMaxDim> xRaw;
  for (int d = 0; d < spatial_dim_; ++d) {
    xRaw[d].reserve(samples.size());
  }

  for (const auto& s : samples) {
    tRaw.push_back(s.t);
    for (int d = 0; d < spatial_dim_; ++d) {
      xRaw[d].push_back(s.x[d]);
    }
  }

  t_axis_ = unique_sorted(tRaw);
  for (int d = 0; d < spatial_dim_; ++d) {
    x_axis_[d] = unique_sorted(xRaw[d]);
  }

  std::size_t expected = t_axis_.size();
  for (int d = 0; d < spatial_dim_; ++d) expected *= x_axis_[d].size();
  if (expected != samples.size()) {
    throw std::runtime_error("TabulatedFunction: CSV file is not a structured grid (expected " +
                             std::to_string(expected) + ", received " +
                             std::to_string(samples.size()) + ")");
  }

  auto is_uniform = [](const std::vector<double>& axis) {
    if (axis.size() < 3) return true;
    const double h = axis[1] - axis[0];
    const double tol = 1e-9 * std::max(1.0, std::abs(h));
    for (std::size_t i = 2; i < axis.size(); ++i) {
      if (std::abs((axis[i] - axis[i - 1]) - h) > tol) return false;
    }
    return true;
  };
  t_uniform_ = is_uniform(t_axis_);
  for (int d = 0; d < spatial_dim_; ++d) x_uniform_[d] = is_uniform(x_axis_[d]);

  data_.resize(expected);

  auto idx_of = [](const std::vector<double>& axis, double val) {
    return static_cast<int>(std::lower_bound(axis.begin(), axis.end(), val) - axis.begin());
  };

  std::vector<bool> filled(expected, false);

  for (const auto& s : samples) {
    const int it = idx_of(t_axis_, s.t);
    std::array<int, kMaxDim> ix{};
    for (int d = 0; d < spatial_dim_; ++d) ix[d] = idx_of(x_axis_[d], s.x[d]);
    const std::size_t k = flat_index(it, ix);
    if (filled[k]) {
      std::ostringstream oss;
      oss << "TabulatedFunction: duplicate grid point in CSV file (t = " << s.t;
      for (int d = 0; d < spatial_dim_; ++d) oss << ", x" << d + 1 << " = " << s.x[d];
      oss << ")";
      throw std::runtime_error(oss.str());
    }
    filled[k] = true;
    data_[k] = s.F;
  }
}

/**
 * @brief Convert a (time, inputs) multi-index into a flat index into data_.
 *
 * The data is stored in row-major order with time as the slowest-varying
 * index and the last input dimension as the fastest:
 *
 *   flat = ((it * N0 + i0) * N1 + i1) * N2 + i2   (3 inputs)
 *
 * where Nd = x_axis_[d].size().
 *
 * @param it  Time index, in [0, t_axis_.size() - 1].
 * @param idx Input indices; only the first spatial_dim_ entries are used,
 *            each idx[d] in [0, x_axis_[d].size() - 1].
 * @return std::size_t Position of the corresponding value in data_.
 *
 * @pre it and each idx[d] are within their axis bounds (checked with
 *      MFEM_ASSERT in debug builds only).
 */
std::size_t TabulatedFunction::flat_index(int it, const std::array<int, kMaxDim>& idx) const {
  MFEM_ASSERT(it >= 0 && static_cast<std::size_t>(it) < t_axis_.size(),
              "TabulatedFunction::flat_index: time index out of range");

  std::size_t flat = static_cast<std::size_t>(it);
  for (int d = 0; d < spatial_dim_; ++d) {
    MFEM_ASSERT(idx[d] >= 0 && static_cast<std::size_t>(idx[d]) < x_axis_[d].size(),
                "TabulatedFunction::flat_index: input index out of range");
    flat = flat * x_axis_[d].size() + static_cast<std::size_t>(idx[d]);
  }

  MFEM_ASSERT(flat < data_.size(), "TabulatedFunction::flat_index: flat index exceeds data size");
  return flat;
}

/**
 * @brief Locate a value on a 1D axis and compute its local interpolation weight.
 *
 * Finds the cell [axis[i], axis[i+1]] containing val and returns the lower
 * index together with the normalized position of val within that cell.
 *
 *  - Uniform axis: the index is computed directly as
 *    floor((val - axis[0]) / h), with h = axis[1] - axis[0] (O(1)).
 *  - Non-uniform axis: the index is found by binary search with
 *    std::upper_bound (O(log n)).
 *
 * Out-of-range values are clamped: the index is restricted to [0, n-2] and the
 * weight to [0, 1], so points outside the axis take the boundary value
 * (constant extrapolation).
 *
 * @param axis    Sorted (strictly increasing) grid coordinates along one axis.
 * @param uniform True if the axis is regularly spaced, enabling the O(1) lookup.
 * @param val     Coordinate to locate.
 * @return std::pair<int, double> {i, a} where i is the lower cell index and
 *         a in [0, 1] is the local weight, such that
 *         val ≈ (1 - a) * axis[i] + a * axis[i+1].
 *         For a single-point axis (n == 1), returns {0, 0.0}.
 */
std::pair<int, double> TabulatedFunction::locate(const std::vector<double>& axis, bool uniform,
                                                 double val) const {
  const int n = static_cast<int>(axis.size());
  if (n == 1) return {0, 0.0};

  if (uniform) {
    const double h = axis[1] - axis[0];
    const int i = std::clamp(static_cast<int>(std::floor((val - axis[0]) / h)), 0, n - 2);
    return {i, std::clamp((val - axis[i]) / h, 0.0, 1.0)};
  }
  int i = static_cast<int>(std::upper_bound(axis.begin(), axis.end(), val) - axis.begin()) - 1;
  i = std::clamp(i, 0, n - 2);
  return {i, std::clamp((val - axis[i]) / (axis[i + 1] - axis[i]), 0.0, 1.0)};
}

/**
 * @brief Interpolate the table at a point given in CSV coordinates.
 *
 * The point is expressed in the coordinates of the CSV input columns
 * (x1..xN) as stored in the file
 *
 * Performs a separable multilinear interpolation on the rectilinear
 * (time × inputs) grid:
 *  - the time axis and each input axis are searched independently with
 *    locate(), which returns the lower cell index and the normalized local
 *    coordinate (weight) within that cell;
 *  - at each of the two bracketing time slices, the table is interpolated over
 *    the inputs by summing the contributions of the 2^N cell corners, each
 *    weighted by the product of the per-axis linear weights;
 *  - the two values are then linearly interpolated in time.
 *
 * Degenerate axes are handled explicitly: a single-point input axis has no
 * upper neighbour (its upper corners are skipped), and a single time value
 * gives a stationary field (no time interpolation).
 *
 * @param t   Time at which the table is evaluated.
 * @param pos Point in CSV input coordinates. Its size must be equal to the
 *            number of input columns of the CSV.
 * @return double Interpolated value at (t, pos).
 */
double TabulatedFunction::interpolate(double t, std::span<const double> pos) const {
  MFEM_VERIFY(static_cast<int>(pos.size()) == spatial_dim_,
              "TabulatedFunction: point dimension not compatible with the current table");

  auto [it, at] = locate(t_axis_, t_uniform_, t);

  std::array<int, kMaxDim> ix{};
  std::array<double, kMaxDim> ax{};
  for (int d = 0; d < spatial_dim_; ++d) {
    auto [i, a] = locate(x_axis_[d], x_uniform_[d], pos[d]);
    ix[d] = i;
    ax[d] = a;
  }

  auto local_lerp = [](double a, double b, double w) { return a + w * (b - a); };

  auto spatialAt = [&](int itc) {
    const int nCorners = 1 << spatial_dim_;
    double acc = 0.0;
    for (int c = 0; c < nCorners; ++c) {
      double w = 1.0;
      bool valid = true;
      std::array<int, kMaxDim> corner = ix;
      for (int d = 0; d < spatial_dim_; ++d) {
        const bool hi = (c >> d) & 1;
        // No upper neighbour on a single-point axis
        if (hi && x_axis_[d].size() == 1) {
          valid = false;
          break;
        }
        corner[d] += hi ? 1 : 0;
        w *= hi ? ax[d] : (1.0 - ax[d]);
      }
      if (valid) acc += w * data_[flat_index(itc, corner)];
    }
    return acc;
  };

  // Stationary field (single time value): no time interpolation
  if (t_axis_.size() == 1) return spatialAt(0);
  return local_lerp(spatialAt(it), spatialAt(it + 1), at);
}

/**
 * @brief Evaluate the table at a point in CSV input coordinates given as an mfem::Vector.
 *
 * Convenience adapter: wraps the vector data in a non-owning std::span and
 * forwards to interpolate(double, std::span<const double>).
 *
 * @param t   Time at which the table is evaluated.
 * @param pos Point in CSV input coordinates
 * @return double Interpolated value at (t, pos).
 */
double TabulatedFunction::interpolate(double t, const mfem::Vector& pos) const {
  return interpolate(t, std::span<const double>(pos.GetData(), pos.Size()));
}