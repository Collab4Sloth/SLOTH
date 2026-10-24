/**
 * @file test_FunctionCoefficients.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for FunctionCoefficient and its derived classes (W, Fw, H, Log, GradientEnergy, Constant, Sum, Product)
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
#include <algorithm>
#include <cmath>
#include <vector>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Coefficients/CommonCoefficients.hpp"
#include "Coefficients/ConstantCoefficient.hpp"
#include "Coefficients/ProductFunctionCoefficients.hpp"
#include "Coefficients/SumFunctionCoefficients.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
constexpr double tol = 1.e-12;

/**
 * @brief Checks the gradient and the Hessian of a coefficient against centred
 * finite differences, at the point x:
 *   dF/dx_i       ~ (F(x + h e_i) - F(x - h e_i)) / 2h
 *   d2F/dx_i dx_j ~ (dF/dx_i(x + h e_j) - dF/dx_i(x - h e_j)) / 2h
 *
 * The finite-difference error is O(h^2) (truncation) + O(eps/h) (round-off):
 * with h = 1e-5, a tolerance of 1e-6 (relative to max(1, |value|)) is safe,
 * while a wrong derivative is off by O(1).
 *
 * The dimension is always given: without it, FunctionCoefficient passes -1,
 * converted to a huge unsigned int (see the note at the end of this file).
 */
void check_derivatives(FunctionCoefficient& c, const std::vector<double>& x, int dim = 1) {
  const double h = 1.e-5;
  const auto n = static_cast<int>(x.size());
  auto shifted = [&](int k, double dx) {
    std::vector<double> y = x;
    y[static_cast<size_t>(k)] += dx;
    return y;
  };
  // Matcher: |value - expected| <= 1e-6 * max(1, |expected|)
  auto near = [](double expected) {
    return WithinAbs(expected, 1.e-6 * std::max(1., std::abs(expected)));
  };

  for (int i = 0; i < n; ++i) {
    CAPTURE(i);
    const double fd = (c.eval_f(shifted(i, h), dim) - c.eval_f(shifted(i, -h), dim)) / (2. * h);
    CHECK_THAT(c.eval_gradient(i, x, dim), near(fd));
    for (int j = 0; j < n; ++j) {
      CAPTURE(j);
      const double fd2 = (c.eval_gradient(i, shifted(j, h), dim) -
                          c.eval_gradient(i, shifted(j, -h), dim)) /
                         (2. * h);
      CHECK_THAT(c.eval_hessian(i, j, x, dim), near(fd2));
    }
  }
}
}  // namespace

// -----------------------------------------------------------------------------
// Derivatives of every coefficient: one generic check, run for each type
// (TEMPLATE_TEST_CASE with types) and several points (GENERATE).
// A new coefficient is checked by adding its type to the list.
// -----------------------------------------------------------------------------
TEMPLATE_TEST_CASE("FunctionCoefficient: gradient and Hessian match finite differences",
                   "[Coefficients][derivatives]", W, Fw, H, Log) {
  // GENERATE: the test case is run once per value
  const double x = GENERATE(0.1, 0.3, 0.5, 0.8);
  CAPTURE(x);
  TestType c(2.5);  // non-unit prefactor, so that a forgotten prefactor is detected
  check_derivatives(c, {x});
}

// -----------------------------------------------------------------------------
// Values: formulas written from the documentation of each potential
// -----------------------------------------------------------------------------
TEST_CASE("CommonCoefficients: W = x^2 (1-x)^2 (double well)", "[Coefficients]") {
  W w;  // default prefactor: 1
  CHECK_THAT(w.eval_f(std::vector{0.}, 1), WithinAbs(0., tol));  // minima in 0 and 1
  CHECK_THAT(w.eval_f(std::vector{1.}, 1), WithinAbs(0., tol));
  CHECK_THAT(w.eval_f(std::vector{0.5}, 1), WithinRel(0.0625, tol));  // barrier
  CHECK_THAT(w.eval_gradient(0, std::vector{0.5}, 1), WithinAbs(0., tol));

  W w3(3.);
  CHECK_THAT(w3.eval_f(std::vector{0.5}, 1), WithinRel(3. * 0.0625, tol));
}

TEST_CASE("CommonCoefficients: Fw = (x^2 - 1)^2 / 4 (double well in -1, 1)", "[Coefficients]") {
  Fw fw;
  CHECK_THAT(fw.eval_f(std::vector{-1.}, 1), WithinAbs(0., tol));
  CHECK_THAT(fw.eval_f(std::vector{1.}, 1), WithinAbs(0., tol));
  CHECK_THAT(fw.eval_f(std::vector{0.}, 1), WithinRel(0.25, tol));
}

TEST_CASE("CommonCoefficients: H = x^3 (6x^2 - 15x + 10) (interpolation)", "[Coefficients]") {
  H h;
  CHECK_THAT(h.eval_f(std::vector{0.}, 1), WithinAbs(0., tol));
  CHECK_THAT(h.eval_f(std::vector{1.}, 1), WithinRel(1., tol));
  CHECK_THAT(h.eval_f(std::vector{0.5}, 1), WithinRel(0.5, tol));  // symmetric
  // Flat at both ends: H'(0) = H'(1) = 0
  CHECK_THAT(h.eval_gradient(0, std::vector{0.}, 1), WithinAbs(0., tol));
  CHECK_THAT(h.eval_gradient(0, std::vector{1.}, 1), WithinAbs(0., tol));
}

TEST_CASE("CommonCoefficients: Log = x ln x + (1-x) ln(1-x)", "[Coefficients]") {
  Log l;
  CHECK_THAT(l.eval_f(std::vector{0.5}, 1), WithinRel(std::log(0.5), tol));
  CHECK_THAT(l.eval_gradient(0, std::vector{0.5}, 1), WithinAbs(0., tol));
  CHECK_THAT(l.eval_hessian(0, 0, std::vector{0.5}, 1), WithinRel(4., tol));

  SECTION("x is clamped to [1e-10, 1 - 1e-10]: finite values at 0 and 1") {
    CHECK(std::isfinite(l.eval_f(std::vector{0.}, 1)));
    CHECK(std::isfinite(l.eval_gradient(0, std::vector{1.}, 1)));
    CHECK(std::isfinite(l.eval_hessian(0, 0, std::vector{-0.5}, 1)));
  }
}

TEST_CASE("CommonCoefficients: GradientEnergy = |g|^2 / 2", "[Coefficients]") {
  GradientEnergy ge(2.);
  // Input: the DIM components of the gradient
  CHECK_THAT(ge.eval_f(std::vector{3., 4.}, 2), WithinRel(2. * 12.5, tol));
  CHECK_THAT(ge.eval_f(std::vector{3., 4., 12.}, 3), WithinRel(2. * 84.5, tol));
  // Only the first DIM components are used
  CHECK_THAT(ge.eval_f(std::vector{3., 4., 12.}, 1), WithinRel(2. * 4.5, tol));
  // Gradient and Hessian with respect to the variable (not to its gradient): 0
  CHECK_THAT(ge.eval_gradient(0, std::vector{3., 4.}, 2), WithinAbs(0., tol));
  CHECK_THAT(ge.eval_hessian(0, 0, std::vector{3., 4.}, 2), WithinAbs(0., tol));
}

TEST_CASE("ConstantCoefficient", "[Coefficients]") {
  ConstantCoefficient c(7.5);
  CHECK_THAT(c.eval_f(std::vector{0.3}, 1), WithinRel(7.5, tol));
  CHECK_THAT(c.eval_f(std::vector{123.}, 1), WithinRel(7.5, tol));  // input ignored
  CHECK_THAT(c.eval_gradient(0, std::vector{0.3}, 1), WithinAbs(0., tol));
  CHECK_THAT(c.eval_hessian(0, 0, std::vector{0.3}, 1), WithinAbs(0., tol));

  SECTION("several variables: n x n Hessian") {
    const std::vector<double> x{0.3, 0.4, 0.5};
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        CHECK(c.eval_hessian(i, j, x, 1) == 0.);
      }
    }
    check_derivatives(c, x);
  }
}

// -----------------------------------------------------------------------------
// FunctionCoefficient: common interface
// -----------------------------------------------------------------------------
TEST_CASE("FunctionCoefficient: evaluation overloads", "[Coefficients]") {
  // The auxiliary and explicit values are ignored by W: the three overloads
  // must give the same result.
  W w(2.);
  const std::vector<double> x{0.3};
  const std::vector<double> exp_x{0.9};
  const std::vector<double> aux{1000.};

  const double f = w.eval_f(x, 1);
  CHECK_THAT(w.eval_f(x, aux, 1), WithinRel(f, tol));
  CHECK_THAT(w.eval_f(x, exp_x, aux, 1), WithinRel(f, tol));

  const double g = w.eval_gradient(0, x, 1);
  CHECK_THAT(w.eval_gradient(0, x, aux, 1), WithinRel(g, tol));
  CHECK_THAT(w.eval_gradient(0, x, exp_x, aux, 1), WithinRel(g, tol));

  const double hs = w.eval_hessian(0, 0, x, 1);
  CHECK_THAT(w.eval_hessian(0, 0, x, aux, 1), WithinRel(hs, tol));
  CHECK_THAT(w.eval_hessian(0, 0, x, exp_x, aux, 1), WithinRel(hs, tol));
}

TEST_CASE("FunctionCoefficient: time and coordinates", "[Coefficients]") {
  W w;
  CHECK(w.get_time() == 0.);
  w.set_time(3.5);
  CHECK(w.get_time() == 3.5);
  CHECK_FALSE(w.uses_coordinates());
  CHECK_NOTHROW(w.set_coordinates(std::vector{0.1, 0.2}));
}

// -----------------------------------------------------------------------------
// Sum and product of coefficients (non-owning: the terms must outlive them)
// -----------------------------------------------------------------------------
TEST_CASE("SumCoefficient", "[Coefficients]") {
  W w(2.);
  H h(3.);
  SumCoefficient sum{&w, &h};

  const std::vector<double> x{0.3};
  CHECK_THAT(sum.eval_f(x, 1), WithinRel(w.eval_f(x, 1) + h.eval_f(x, 1), tol));
  check_derivatives(sum, x);

  SECTION("set_time is propagated to every term") {
    sum.set_time(2.);
    CHECK(sum.get_time() == 2.);
    CHECK(w.get_time() == 2.);
    CHECK(h.get_time() == 2.);
  }
  SECTION("coordinates") {
    CHECK_FALSE(sum.uses_coordinates());
    CHECK_NOTHROW(sum.set_coordinates(std::vector{0.1, 0.2}));
  }
}

TEST_CASE("ProductCoefficient", "[Coefficients]") {
  W w(2.);
  Log l(3.);
  ProductCoefficient prod{&w, &l};

  const std::vector<double> x{0.3};
  CHECK_THAT(prod.eval_f(x, 1), WithinRel(w.eval_f(x, 1) * l.eval_f(x, 1), tol));
  check_derivatives(prod, x);

  SECTION("set_time is propagated to every term") {
    prod.set_time(2.);
    CHECK(w.get_time() == 2.);
    CHECK(l.get_time() == 2.);
  }
  SECTION("coordinates") {
    CHECK_FALSE(prod.uses_coordinates());
    CHECK_NOTHROW(prod.set_coordinates(std::vector{0.1, 0.2}));
  }
}

// Derivatives by the product rule: defined where a factor vanishes.
// W and H both vanish at x = 0, with W'(0) = H'(0) = 0:
//   (W H)'(0) = 0,  (W H)''(0) = W''(0) H(0) + 2 W'(0) H'(0) + W(0) H''(0) = 0
// W(1) = W'(1) = 0: (W H)'(1) = 0, (W H)''(1) = W''(1) H(1).
TEST_CASE("ProductCoefficient: derivatives where a factor vanishes", "[Coefficients]") {
  W w;
  H h;
  ProductCoefficient prod{&w, &h};
  const std::vector<double> zero{0.};
  CHECK_THAT(prod.eval_gradient(0, zero, 1), WithinAbs(0., tol));
  CHECK_THAT(prod.eval_hessian(0, 0, zero, 1), WithinAbs(0., tol));

  SECTION("x = 1") {
    const std::vector<double> one{1.};
    CHECK_THAT(prod.eval_gradient(0, one, 1), WithinAbs(0., tol));
    CHECK_THAT(prod.eval_hessian(0, 0, one, 1),
               WithinAbs(w.eval_hessian(0, 0, one, 1) * h.eval_f(one, 1), tol));
  }
  SECTION("three factors, with and without a vanishing factor") {
    ProductCoefficient prod3{&w, &h, &h};
    check_derivatives(prod3, std::vector{0.5});
    check_derivatives(prod3, std::vector{0.});
  }
}

// Note on the dimension argument: eval_f(x) without dimension passes -1 to
// F(), whose parameter is an unsigned int: it becomes 4294967295. Harmless for
// coefficients that ignore it, but GradientEnergy would then read far beyond
// the input. Not tested here (undefined behaviour).
