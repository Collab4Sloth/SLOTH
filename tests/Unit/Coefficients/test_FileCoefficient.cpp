/**
 * @file test_FileCoefficient.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for FileCoefficient (coefficient read from a CSV file)
 * @version 0.1
 * @date 2026-10-22
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
#include <array>
#include <cmath>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Coefficients/Coefficient.hpp"
#include "Coefficients/FileCoefficient.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
constexpr double tol = 1.e-10;

/**
 * @brief CSV content "t,<names>,F" on the grid t_axis x axes[0] x ... ,
 * with F computed by f(t, x).
 */
std::string make_csv(const std::string& names, const std::vector<double>& t_axis,
                     const std::vector<std::vector<double>>& axes,
                     const std::function<double(double, const std::vector<double>&)>& f) {
  std::ostringstream out;
  out.precision(17);
  out << "t," << names << ",F\n";

  // Number of spatial nodes, then node k -> indices on each axis
  // (the last axis varies fastest, like nested loops)
  std::size_t nb_nodes = 1;
  for (const auto& axis : axes) nb_nodes *= axis.size();

  for (const double t : t_axis) {
    for (std::size_t k = 0; k < nb_nodes; ++k) {
      std::vector<double> x(axes.size());
      std::size_t rest = k;
      for (std::size_t d = axes.size(); d-- > 0;) {
        x[d] = axes[d][rest % axes[d].size()];
        rest /= axes[d].size();
      }
      out << t;
      for (const double xd : x) out << "," << xd;
      out << "," << f(t, x) << "\n";
    }
  }
  return out.str();
}

/// @brief Evaluates a Coefficient at time t and position x
double eval(Coefficient& c, double t, const std::vector<double>& x) {
  c.set_time(t);
  return c.compute(x, std::span<const double>{}, std::span<const double>{},
                   static_cast<int>(x.size()));
}

/// @brief Reference affine field: F = t + 2x + 3y + 4z
double affine3d(double t, const std::vector<double>& x) {
  return t + 2. * x[0] + 3. * x[1] + 4. * x[2];
}
}  // namespace

// -----------------------------------------------------------------------------
// Affine field on a 3D + time grid.
// Multilinear interpolation is exact for an affine field, at any point.
// -----------------------------------------------------------------------------
TEST_CASE("FileCoefficient: 3D + time, affine field", "[Coefficients][file]") {
  const sloth_unit::TemporaryFile csv(
      "affine3d.csv", make_csv("x,y,z", {0., 1., 2.},
                               {{0., 1., 2., 3.}, {0., 1., 2., 3.}, {0., 1., 2., 3.}}, affine3d));
  Coefficient c{FileCoefficient(csv.path())};

  CHECK(c.is_file());
  CHECK(c.get_type() == GlossaryType::System);  // Glossary::Default without quantity

  SECTION("grid node") {
    CHECK_THAT(eval(c, 1., {2., 1., 3.}), WithinRel(affine3d(1., {2., 1., 3.}), tol));
  }
  SECTION("inside a cell, between two times") {
    CHECK_THAT(eval(c, 1.5, {1.3, 2.7, 0.4}), WithinRel(affine3d(1.5, {1.3, 2.7, 0.4}), tol));
  }
  SECTION("on the boundary of the grid") {
    CHECK_THAT(eval(c, 2., {3., 3., 3.}), WithinRel(affine3d(2., {3., 3., 3.}), tol));
  }
  SECTION("outside the grid: constant extrapolation (clamped to the boundary)") {
    CHECK_THAT(eval(c, 5., {3., 3., 10.}), WithinRel(affine3d(2., {3., 3., 3.}), tol));
    CHECK_THAT(eval(c, -1., {-1., 0., 0.}), WithinAbs(affine3d(0., {0., 0., 0.}), tol));
  }
}

TEST_CASE("FileCoefficient: 1D and 2D fields", "[Coefficients][file]") {
  SECTION("1D") {
    auto f = [](double t, const std::vector<double>& x) { return 10. * t - x[0]; };
    const sloth_unit::TemporaryFile csv("f1d.csv", make_csv("x", {0., 1.}, {{0., 0.5, 1.}}, f));
    FileCoefficient fc(csv.path());
    CHECK(fc.input_dim() == 1);
    CHECK_THAT(fc.value(0.5, std::vector{0.25}), WithinRel(f(0.5, {0.25}), tol));
  }
  SECTION("2D, non-uniform grid") {
    auto f = [](double t, const std::vector<double>& x) { return t + x[0] - 2. * x[1]; };
    const sloth_unit::TemporaryFile csv(
        "f2d.csv", make_csv("x,y", {0., 1.}, {{0., 0.1, 1.}, {0., 0.7, 0.8, 2.}}, f));
    FileCoefficient fc(csv.path());
    CHECK(fc.input_dim() == 2);
    CHECK_THAT(fc.value(0.3, std::vector{0.4, 1.5}), WithinRel(f(0.3, {0.4, 1.5}), tol));
  }
}

TEST_CASE("FileCoefficient: Coefficient with a physical quantity", "[Coefficients][file]") {
  const sloth_unit::TemporaryFile csv(
      "temperature.csv",
      make_csv("x", {0., 1.}, {{0., 1.}},
               [](double t, const std::vector<double>& x) { return 300. + 100. * t + x[0]; }));
  Coefficient c(Glossary::Temperature, FileCoefficient(csv.path()));
  CHECK(c.is_file());
  CHECK(c.get_type() == GlossaryType::Temperature);
  CHECK_THAT(eval(c, 0.5, {0.5}), WithinRel(350.5, tol));
}

TEST_CASE("FileCoefficient: axisymmetric geometry", "[Coefficients][file]") {
  // The CSV is given in (r) or (r, z); the coefficient is evaluated at mesh
  // coordinates (x, y[, z]) with r = |(x, y) - axis_center|.
  SECTION("(t, r, F), evaluated in 2D") {
    auto f = [](double t, const std::vector<double>& r) { return t + 2. * r[0]; };
    const sloth_unit::TemporaryFile csv("axi_r.csv",
                                        make_csv("r", {0., 1.}, {{0., 5., 10.}}, f));
    FileCoefficient fc(csv.path(), Geometry::Axisymmetric, {1., 1.});
    CHECK(fc.geometry() == Geometry::Axisymmetric);
    CHECK(fc.axis_center() == std::array{1., 1.});
    // (4, 5) is at distance 5 from (1, 1)
    CHECK_THAT(fc.value(0.5, std::vector{4., 5.}), WithinRel(f(0.5, {5.}), tol));
  }
  SECTION("(t, r, z, F), evaluated in 3D") {
    auto f = [](double t, const std::vector<double>& rz) { return t + rz[0] + 3. * rz[1]; };
    const sloth_unit::TemporaryFile csv(
        "axi_rz.csv", make_csv("r,z", {0., 1.}, {{0., 2., 4.}, {0., 1.}}, f));
    FileCoefficient fc(csv.path(), Geometry::Axisymmetric);  // axis through (0, 0)
    CHECK_THAT(fc.value(1., std::vector{0., 3., 0.5}), WithinRel(f(1., {3., 0.5}), tol));
  }
  SECTION("errors") {
    const sloth_unit::TemporaryFile csv3d(
        "axi_3d.csv", make_csv("x,y,z", {0., 1.}, {{0., 1.}, {0., 1.}, {0., 1.}}, affine3d));
    SLOTH_CHECK_MFEM_ERROR_WITH(FileCoefficient(csv3d.path(), Geometry::Axisymmetric),
                                "axisymmetric CSV must be");

    const sloth_unit::TemporaryFile csv_rz(
        "axi_rz2.csv",
        make_csv("r,z", {0., 1.}, {{0., 1.}, {0., 1.}},
                 [](double, const std::vector<double>&) { return 0.; }));
    FileCoefficient fc(csv_rz.path(), Geometry::Axisymmetric);
    // (r, z) field needs 3 mesh coordinates
    SLOTH_CHECK_MFEM_ERROR_WITH(fc.value(0., std::vector{1., 1.}), "mesh dimension too small");
  }
}

TEST_CASE("FileCoefficient: errors", "[Coefficients][file]") {
  SECTION("missing file") {
    SLOTH_CHECK_MFEM_ERROR_WITH(FileCoefficient("/nonexistent/sloth_unit.csv"),
                                "Error while trying to open");
  }
  SECTION("no gradient nor Hessian for a tabulated field") {
    const sloth_unit::TemporaryFile csv(
        "nograd.csv", make_csv("x", {0., 1.}, {{0., 1.}},
                               [](double, const std::vector<double>& x) { return x[0]; }));
    FileCoefficient fc(csv.path());
    SLOTH_CHECK_MFEM_ERROR_WITH(fc.eval_gradient(0, std::vector{0.5}, 1), "gradient");
    SLOTH_CHECK_MFEM_ERROR_WITH(fc.eval_hessian(0, 0, std::vector{0.5}, 1), "Hessian");
  }
}
