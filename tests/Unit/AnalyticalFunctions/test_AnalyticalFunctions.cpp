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

TEMPLATE_TEST_CASE_SIG("AnalyticalFunctions: Uniform", "[AnalyticalFunctions]", ((int DIM), DIM),
                       1, 2, 3) {
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
  const AnalyticalFunctions<1> af(AnalyticalFunctionsType::HyperbolicTangent, center, a,
                                  thickness, radius);
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
