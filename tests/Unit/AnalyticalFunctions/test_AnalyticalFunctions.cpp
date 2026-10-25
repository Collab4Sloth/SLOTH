/**
 * @file test_AnalyticalFunctions.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for AnalyticalFunctions<DIM>
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
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <functional>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "Options/AnalyticalFunctionOptions.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
/// @brief Point of dimension DIM: (x0, 0, ..., 0)
template <int DIM>
mfem::Vector point_x(double x0) {
  mfem::Vector x(DIM);
  x = 0.;
  x[0] = x0;
  return x;
}
}  // namespace

// -----------------------------------------------------------------------------
// One test, written once, run for DIM = 1, 2 and 3.
// TEMPLATE_TEST_CASE_SIG is the Catch2 version for non-type template
// parameters: ((int DIM), DIM) declares the parameter, then 1, 2, 3 are the
// values. Catch2 creates three test cases: "... - 1", "... - 2", "... - 3".
// -----------------------------------------------------------------------------
TEMPLATE_TEST_CASE_SIG("AnalyticalFunctions: Heaviside, step at x = 0.5", "[AnalyticalFunctions]",
                       ((int DIM), DIM), 1, 2, 3) {
  const AnalyticalFunctions<DIM> af(AnalyticalFunctionsType::Heaviside);
  const auto f = af.getFunction();

  CHECK(f(point_x<DIM>(0.25), 0.) == 0.);
  CHECK(f(point_x<DIM>(0.5), 0.) == 1.);  // 1 at the step itself
  CHECK(f(point_x<DIM>(0.75), 0.) == 1.);
  CHECK(f(point_x<DIM>(0.75), 10.) == 1.);  // time-independent
}

TEMPLATE_TEST_CASE_SIG("AnalyticalFunctions: Uniform", "[AnalyticalFunctions]", ((int DIM), DIM), 1,
                       2, 3) {
  const AnalyticalFunctions<DIM> af(AnalyticalFunctionsType::Uniform, 273.15);
  const auto f = af.getFunction();

  CHECK_THAT(f(point_x<DIM>(0.), 0.), WithinRel(273.15, 1.e-12));
  CHECK_THAT(f(point_x<DIM>(3.), 5.), WithinRel(273.15, 1.e-12));

  SECTION("exactly one argument is expected") {
    SLOTH_CHECK_MFEM_ERROR_WITH(AnalyticalFunctions<DIM>(AnalyticalFunctionsType::Uniform, 1., 2.),
                                "only one argument");
  }
}

// -----------------------------------------------------------------------------
// Dimension-specific behaviour: ordinary test cases
// -----------------------------------------------------------------------------
TEST_CASE("AnalyticalFunctions: HyperbolicTangent 1D", "[AnalyticalFunctions]") {
  // Arguments: center_x, a_x, thickness, radius
  const double center = 0.5, a = 1., thickness = 0.01, radius = 0.2;
  const AnalyticalFunctions<1> af(AnalyticalFunctionsType::HyperbolicTangent, center, a, thickness,
                                  radius);
  const auto f = af.getFunction();

  // 0.5 on the interface, 0 inside, 1 outside
  CHECK_THAT(f(point_x<1>(center + radius), 0.), WithinRel(0.5, 1.e-12));
  CHECK_THAT(f(point_x<1>(center), 0.), WithinAbs(0., 1.e-10));
  CHECK_THAT(f(point_x<1>(center + 2. * radius), 0.), WithinAbs(1., 1.e-10));

  SECTION("four arguments are expected") {
    SLOTH_CHECK_MFEM_ERROR_WITH(
        AnalyticalFunctions<1>(AnalyticalFunctionsType::HyperbolicTangent, 0.5, 1., 0.01),
        "four arguments");
  }
}

TEST_CASE("AnalyticalFunctionsType::from", "[AnalyticalFunctions]") {
  CHECK(AnalyticalFunctionsType::from("Heaviside") == AnalyticalFunctionsType::Heaviside);
  CHECK(AnalyticalFunctionsType::from("uniform") == AnalyticalFunctionsType::Uniform);
  CHECK(AnalyticalFunctionsType::from("PARABOLIC") == AnalyticalFunctionsType::Parabolic);

  SECTION("unknown keyword") {
    SLOTH_CHECK_MFEM_ERROR_WITH(AnalyticalFunctionsType::from("Gaussian"), "not a valid value");
  }
}

// -----------------------------------------------------------------------------
// HyperbolicTangent in 2D and 3D: circular / spherical interface
// -----------------------------------------------------------------------------
TEST_CASE("AnalyticalFunctions: HyperbolicTangent 2D", "[AnalyticalFunctions]") {
  // Arguments: center_x, center_y, a_x, a_y, thickness, radius
  const AnalyticalFunctions<2> af(AnalyticalFunctionsType::HyperbolicTangent, 0.5, 0.5, 1., 1.,
                                  0.01, 0.25);
  const auto f = af.getFunction();

  // 0.5 on the circle of radius 0.25 centred in (0.5, 0.5), whatever the direction
  const double c = 0.25 / std::sqrt(2.);
  CHECK_THAT(f(sloth_unit::point({0.75, 0.5}), 0.), WithinRel(0.5, 1.e-12));
  CHECK_THAT(f(sloth_unit::point({0.5 + c, 0.5 + c}), 0.), WithinRel(0.5, 1.e-12));
  CHECK_THAT(f(sloth_unit::point({0.5, 0.5}), 0.), WithinAbs(0., 1.e-10));  // centre
  CHECK_THAT(f(sloth_unit::point({1., 1.}), 0.), WithinAbs(1., 1.e-10));    // far outside

  SECTION("six arguments are expected") {
    SLOTH_CHECK_MFEM_ERROR_WITH(
        AnalyticalFunctions<2>(AnalyticalFunctionsType::HyperbolicTangent, 0.5, 0.5, 1., 1.),
        "getHyperbolicTangent");
  }
}

TEST_CASE("AnalyticalFunctions: HyperbolicTangent 3D", "[AnalyticalFunctions]") {
  // Arguments: center_x, center_y, center_z, a_x, a_y, a_z, thickness, radius
  const AnalyticalFunctions<3> af(AnalyticalFunctionsType::HyperbolicTangent, 0., 0., 0., 1., 1.,
                                  1., 0.01, 1.);
  const auto f = af.getFunction();

  // 0.5 on the unit sphere
  const double c = 1. / std::sqrt(3.);
  CHECK_THAT(f(sloth_unit::point({0., 0., 1.}), 0.), WithinRel(0.5, 1.e-12));
  CHECK_THAT(f(sloth_unit::point({c, c, c}), 0.), WithinRel(0.5, 1.e-12));
  CHECK_THAT(f(sloth_unit::point({0., 0., 0.}), 0.), WithinAbs(0., 1.e-10));

  SECTION("eight arguments are expected") {
    SLOTH_CHECK_MFEM_ERROR_WITH(
        AnalyticalFunctions<3>(AnalyticalFunctionsType::HyperbolicTangent, 0., 0., 0., 1.),
        "getHyperbolicTangent");
  }
}

// -----------------------------------------------------------------------------
// Sinusoide and Sinusoide2 (2D only): manufactured solutions
// -----------------------------------------------------------------------------
TEST_CASE("AnalyticalFunctions: Sinusoide 2D", "[AnalyticalFunctions]") {
  const double k = 3.;
  const auto x = sloth_unit::point({0.3, 0.4});
  const double t = 0.7;
  // u = exp(-2t) sin(x + y)
  const double u = std::exp(-2. * t) * std::sin(0.3 + 0.4);

  SECTION("Sinusoide = k u") {
    const AnalyticalFunctions<2> af(AnalyticalFunctionsType::Sinusoide, k);
    CHECK_THAT(af.getFunction()(x, t), WithinRel(k * u, 1.e-12));
  }
  SECTION("Sinusoide2 = k (u^3 - u)") {
    const AnalyticalFunctions<2> af(AnalyticalFunctionsType::Sinusoide2, k);
    CHECK_THAT(af.getFunction()(x, t), WithinRel(k * (u * u * u - u), 1.e-12));
  }
  SECTION("one argument is expected") {
    SLOTH_CHECK_MFEM_ERROR_WITH(AnalyticalFunctions<2>(AnalyticalFunctionsType::Sinusoide),
                                "getSinusoide");
    SLOTH_CHECK_MFEM_ERROR_WITH(AnalyticalFunctions<2>(AnalyticalFunctionsType::Sinusoide2, 1., 2.),
                                "getSinusoide2");
  }
}

// -----------------------------------------------------------------------------
// Parabolic (1D and 2D): temperature profile in a pellet
//   T(r) = f0 + P (rmax^2 - r^2) / (4 pi lambda rmax^2),  r = x[0]
// -----------------------------------------------------------------------------
TEMPLATE_TEST_CASE_SIG("AnalyticalFunctions: Parabolic", "[AnalyticalFunctions]", ((int DIM), DIM),
                       1, 2) {
  // Arguments: rmax, f0, linear power, conductivity
  const double rmax = 4.e-3, f0 = 600., plin = 2.e4, cond = 3.;
  const AnalyticalFunctions<DIM> af(AnalyticalFunctionsType::Parabolic, rmax, f0, plin, cond);
  const auto f = af.getFunction();

  CHECK_THAT(f(point_x<DIM>(rmax), 0.), WithinRel(f0, 1.e-12));  // surface
  CHECK_THAT(f(point_x<DIM>(0.), 0.), WithinRel(f0 + plin / (4. * M_PI * cond), 1.e-12));  // centre

  SECTION("four arguments are expected") {
    SLOTH_CHECK_MFEM_ERROR_WITH(
        AnalyticalFunctions<DIM>(AnalyticalFunctionsType::Parabolic, rmax, f0, plin),
        "getParabolic: four arguments are expected");
  }
}

// -----------------------------------------------------------------------------
// Functions not available in a given dimension
// -----------------------------------------------------------------------------
TEST_CASE("AnalyticalFunctions: not implemented in this dimension", "[AnalyticalFunctions]") {
  SLOTH_CHECK_MFEM_ERROR_WITH(AnalyticalFunctions<1>(AnalyticalFunctionsType::Sinusoide, 1.),
                              "Not implemented");
  SLOTH_CHECK_MFEM_ERROR_WITH(AnalyticalFunctions<1>(AnalyticalFunctionsType::Sinusoide2, 1.),
                              "Not implemented");
  SLOTH_CHECK_MFEM_ERROR_WITH(AnalyticalFunctions<3>(AnalyticalFunctionsType::Sinusoide, 1.),
                              "Not implemented");
  SLOTH_CHECK_MFEM_ERROR_WITH(AnalyticalFunctions<3>(AnalyticalFunctionsType::Sinusoide2, 1.),
                              "Not implemented");
  SLOTH_CHECK_MFEM_ERROR_WITH(
      AnalyticalFunctions<3>(AnalyticalFunctionsType::Parabolic, 1., 1., 1., 1.),
      "Not implemented");
}

// -----------------------------------------------------------------------------
// User-defined function
// -----------------------------------------------------------------------------
TEST_CASE("AnalyticalFunctions: user-defined function", "[AnalyticalFunctions]") {
  const AnalyticalFunctions<2> af(std::function<double(const mfem::Vector&, double)>(
      [](const mfem::Vector& x, double t) { return x[0] * x[1] + t; }));
  CHECK_THAT(af.getFunction()(sloth_unit::point({2., 3.}), 1.), WithinRel(7., 1.e-12));
}
