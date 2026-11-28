/**
 * @file test_AverageCoefficient.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for QuadRule and AverageCoefficient (average of a coefficient over
 * the segment [u_n, u_{n+1}]: midpoint, trapezoid, average vector field)
 * @version 0.1
 * @date 2026-11-27
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
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <functional>
#include <numeric>
#include <span>
#include <vector>

#include "Coefficients/AverageCoefficient.hpp"
#include "Coefficients/CommonCoefficients.hpp"
#include "sloth_unit.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
constexpr double tol = 1.e-12;

/// Matcher: |value - expected| <= t * max(1, |expected|)
auto near(const double expected, const double t = tol) {
  return WithinAbs(expected, t * std::max(1., std::abs(expected)));
}

/**
 * @brief Checks the derivatives of c with respect to u_{n+1} (values x, explicit values xn
 * fixed) against centred finite differences:
 *   - Value: dF/dx_i ~ (F(x + h e_i) - F(x - h e_i)) / 2h;
 *   - both: d2F/dx_i dx_j ~ (dF/dx_i(x + h e_j) - dF/dx_i(x - h e_j)) / 2h.
 * In Gradient mode, the first derivatives are not the derivatives of the value (g(u_{n+1})):
 * only the second derivatives are checked by finite differences.
 */
void check_derivatives(AverageCoefficient& c, const std::vector<double>& x,
                       const std::vector<double>& xn) {
  const double h = 1.e-5;
  const std::vector<double> aux;
  const auto n = static_cast<int>(x.size());
  auto shifted = [&](int k, double dx) {
    std::vector<double> y = x;
    y[static_cast<std::size_t>(k)] += dx;
    return y;
  };
  for (int i = 0; i < n; ++i) {
    CAPTURE(i);
    if (c.average() == Average::Value) {
      const double fd =
          (c.eval_f(shifted(i, h), xn, aux, n) - c.eval_f(shifted(i, -h), xn, aux, n)) / (2. * h);
      CHECK_THAT(c.eval_gradient(i, x, xn, aux, n), near(fd, 1.e-6));
    }
    for (int j = 0; j < n; ++j) {
      CAPTURE(j);
      const double fd2 = (c.eval_gradient(i, shifted(j, h), xn, aux, n) -
                          c.eval_gradient(i, shifted(j, -h), xn, aux, n)) /
                         (2. * h);
      CHECK_THAT(c.eval_hessian(i, j, x, xn, aux, n), near(fd2, 1.e-6));
    }
  }
}
}  // namespace

// -----------------------------------------------------------------------------
// Quadrature rules on [0, 1]
// -----------------------------------------------------------------------------
TEST_CASE("QuadRule: predefined rules", "[Coefficients][Average]") {
  const auto e = QuadRule::end_point();
  CHECK(e.s == std::vector{1.});
  CHECK(e.w == std::vector{1.});
  const auto m = QuadRule::midpoint();
  CHECK(m.s == std::vector{0.5});
  CHECK(m.w == std::vector{1.});
  const auto t = QuadRule::trapezoid();
  CHECK(t.s == std::vector{0., 1.});
  CHECK(t.w == std::vector{0.5, 0.5});
}

TEST_CASE("QuadRule: Gauss rule exact up to its order", "[Coefficients][Average]") {
  const int order = GENERATE(0, 1, 2, 3, 5, 7);
  CAPTURE(order);
  const auto rule = QuadRule::gauss(order);
  REQUIRE(rule.s.size() == rule.w.size());
  CHECK_THAT(std::accumulate(rule.w.begin(), rule.w.end(), 0.), WithinRel(1., tol));
  // int_0^1 s^k ds = 1 / (k + 1)
  for (int k = 0; k <= order; ++k) {
    CAPTURE(k);
    double integral = 0.;
    for (std::size_t q = 0; q < rule.s.size(); ++q) {
      CHECK(rule.s[q] >= 0.);
      CHECK(rule.s[q] <= 1.);
      integral += rule.w[q] * std::pow(rule.s[q], k);
    }
    CHECK_THAT(integral, WithinRel(1. / (k + 1), tol));
  }
}

// -----------------------------------------------------------------------------
// Construction
// -----------------------------------------------------------------------------
TEST_CASE("AverageCoefficient: invalid arguments", "[Coefficients][Average]") {
  W w;
  SLOTH_CHECK_MFEM_ERROR(AverageCoefficient(nullptr, QuadRule::midpoint(), Average::Value));
  SLOTH_CHECK_MFEM_ERROR(AverageCoefficient(&w, QuadRule{{}, {}}, Average::Value));
  SLOTH_CHECK_MFEM_ERROR(AverageCoefficient(&w, QuadRule{{0.5}, {1., 0.}}, Average::Value));
  SLOTH_CHECK_MFEM_ERROR(AverageCoefficient(&w, QuadRule{{1.5}, {1.}}, Average::Value));
  SLOTH_CHECK_MFEM_ERROR(AverageCoefficient(&w, QuadRule{{0.5}, {0.9}}, Average::Value));
  SLOTH_CHECK_MFEM_ERROR(AverageCoefficient(&w, QuadRule{{0., 1.}, {1.5, -0.5}}, Average::Value));
}

TEST_CASE("AverageCoefficient: values at the beginning of the time step required",
          "[Coefficients][Average]") {
  W w;
  AverageCoefficient c(&w, QuadRule::midpoint(), Average::Value);
  const std::vector<double> x{0.3};
  SLOTH_CHECK_MFEM_ERROR(c.eval_f(x, 1));
  SLOTH_CHECK_MFEM_ERROR(c.eval_gradient(0, x, 1));
  SLOTH_CHECK_MFEM_ERROR(c.eval_hessian(0, 0, x, 1));
}

TEST_CASE("AverageCoefficient: properties forwarded to the averaged coefficient",
          "[Coefficients][Average]") {
  W w;
  AverageCoefficient c(&w, QuadRule::midpoint(), Average::Value);
  CHECK(c.is_single_variable() == w.is_single_variable());
  CHECK(c.uses_coordinates() == w.uses_coordinates());
  CHECK(c.batch_nb_variables() == w.batch_nb_variables());
  c.set_time(2.5);
  CHECK(c.get_time() == 2.5);
  CHECK(w.get_time() == 2.5);
}

// -----------------------------------------------------------------------------
// Average of the value (mobility...)
// -----------------------------------------------------------------------------
TEST_CASE("AverageCoefficient: value at the midpoint", "[Coefficients][Average]") {
  W w(2.);
  AverageCoefficient c(&w, QuadRule::midpoint(), Average::Value);
  const double xn = GENERATE(0.1, 0.4);
  const double x = GENERATE(0.2, 0.9);
  CAPTURE(xn, x);
  const std::vector<double> u{x}, un{xn}, aux;
  const std::vector<double> mid{0.5 * (x + xn)};
  CHECK_THAT(c.eval_f(u, un, aux, 1), near(w.eval_f(mid, 1)));
  CHECK_THAT(c.eval_gradient(0, u, un, aux, 1), near(0.5 * w.eval_gradient(0, mid, 1)));
  CHECK_THAT(c.eval_hessian(0, 0, u, un, aux, 1), near(0.25 * w.eval_hessian(0, 0, mid, 1)));
  check_derivatives(c, u, un);
}

TEST_CASE("AverageCoefficient: trapezoid average of the value", "[Coefficients][Average]") {
  W w(2.);
  AverageCoefficient c(&w, QuadRule::trapezoid(), Average::Value);
  const std::vector<double> u{0.7}, un{0.2}, aux;
  CHECK_THAT(c.eval_f(u, un, aux, 1), near(0.5 * (w.eval_f(u, 1) + w.eval_f(un, 1))));
  // Only the value at u_{n+1} depends on u_{n+1}
  CHECK_THAT(c.eval_gradient(0, u, un, aux, 1), near(0.5 * w.eval_gradient(0, u, 1)));
  CHECK_THAT(c.eval_hessian(0, 0, u, un, aux, 1), near(0.5 * w.eval_hessian(0, 0, u, 1)));
  check_derivatives(c, u, un);
}

TEST_CASE("AverageCoefficient: end point rule gives the implicit value",
          "[Coefficients][Average]") {
  W w(2.);
  const std::vector<double> u{0.7}, un{0.2}, aux;
  for (const auto average : {Average::Value, Average::Gradient}) {
    AverageCoefficient c(&w, QuadRule::end_point(), average);
    CHECK_THAT(c.eval_f(u, un, aux, 1), near(w.eval_f(u, 1)));
    CHECK_THAT(c.eval_gradient(0, u, un, aux, 1), near(w.eval_gradient(0, u, 1)));
    CHECK_THAT(c.eval_hessian(0, 0, u, un, aux, 1), near(w.eval_hessian(0, 0, u, 1)));
  }
}

// -----------------------------------------------------------------------------
// Average of the gradient: average vector field (discrete gradient) of a free energy
// -----------------------------------------------------------------------------
TEST_CASE("AverageCoefficient: average vector field of the double well Fw",
          "[Coefficients][Average]") {
  // Fw = (x^2 - 1)^2 / 4, Fw' = x^3 - x (degree 3: Gauss rule of order 3)
  // AVF: int_0^1 Fw'(a + s (b - a)) ds = (a^3 + a^2 b + a b^2 + b^3) / 4 - (a + b) / 2
  Fw fw;
  AverageCoefficient c(&fw, QuadRule::gauss(3), Average::Gradient);
  const double a = GENERATE(-0.9, -0.2, 0.5);
  const double b = GENERATE(-0.6, 0.1, 1.2);
  CAPTURE(a, b);
  const std::vector<double> u{b}, un{a}, aux;
  const double avf = (a * a * a + a * a * b + a * b * b + b * b * b) / 4. - (a + b) / 2.;
  CHECK_THAT(c.eval_gradient(0, u, un, aux, 1), near(avf));
  // Discrete gradient: avf (b - a) = Fw(b) - Fw(a)
  CHECK_THAT(c.eval_gradient(0, u, un, aux, 1) * (b - a), near(fw.eval_f(u, 1) - fw.eval_f(un, 1)));
  // Value: free energy at u_{n+1}
  CHECK_THAT(c.eval_f(u, un, aux, 1), near(fw.eval_f(u, 1)));
  // d avf / db = (a^2 + 2 a b + 3 b^2) / 4 - 1 / 2
  CHECK_THAT(c.eval_hessian(0, 0, u, un, aux, 1),
             near((a * a + 2. * a * b + 3. * b * b) / 4. - 0.5));
  check_derivatives(c, u, un);
}

TEST_CASE("AverageCoefficient: average vector field when u_{n+1} = u_n",
          "[Coefficients][Average]") {
  // Consistency: implicit derivative, half of the implicit second derivative
  W w(3.);
  AverageCoefficient c(&w, QuadRule::gauss(3), Average::Gradient);
  const std::vector<double> u{0.35}, aux;
  CHECK_THAT(c.eval_gradient(0, u, u, aux, 1), near(w.eval_gradient(0, u, 1)));
  CHECK_THAT(c.eval_hessian(0, 0, u, u, aux, 1), near(0.5 * w.eval_hessian(0, 0, u, 1)));
}

// -----------------------------------------------------------------------------
// Average in time: point s of the segment evaluated at t - (1 - s) dt
// -----------------------------------------------------------------------------
namespace {
/// g(u, t) = t u
class TimeLinear final : public FunctionCoefficient {
 protected:
  std::function<double(const std::span<const double>&, const std::span<const double>&,
                       const std::span<const double>&, const unsigned int)>
  F() override {
    return [this](const std::span<const double>& u, const std::span<const double>&,
                  const std::span<const double>&, const unsigned int) {
      return this->time_ * u[0];
    };
  }
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int)>
  GradientF() override {
    return [this](const std::span<const double>&, const std::span<const double>&,
                  const std::span<const double>&, const unsigned int) {
      return std::vector<double>{this->time_};
    };
  }
  std::function<std::vector<double>(const std::span<const double>&, const std::span<const double>&,
                                    const std::span<const double>&, const unsigned int)>
  HessianF() override {
    return [](const std::span<const double>&, const std::span<const double>&,
              const std::span<const double>&, const unsigned int) {
      return std::vector<double>{0.};
    };
  }
};
}  // namespace

TEST_CASE("AverageCoefficient: average in time over the step", "[Coefficients][Average]") {
  TimeLinear g;
  AverageCoefficient c(&g, QuadRule::gauss(2), Average::Value);
  const double a = 0.4;  // u_n
  const double b = 1.0;  // u_{n+1}
  const std::vector<double> u{b}, un{a}, aux;
  c.set_time(1.);

  SECTION("no time step: evaluation at the time of the coefficient") {
    // int_0^1 1 (a + s (b - a)) ds = (a + b) / 2
    CHECK_THAT(c.eval_f(u, un, aux, 1), near(0.5 * (a + b)));
  }
  SECTION("time step dt = 0.5: t(s) = 0.5 + 0.5 s") {
    c.set_time_step(0.5);
    // int_0^1 (0.5 + 0.5 s) (a + s (b - a)) ds = 3 a / 4 + 5 (b - a) / 12
    CHECK_THAT(c.eval_f(u, un, aux, 1), near(0.75 * a + 5. * (b - a) / 12.));
    // d/db: int_0^1 (0.5 + 0.5 s) s ds = 5 / 12
    CHECK_THAT(c.eval_gradient(0, u, un, aux, 1), near(5. / 12.));
    // Time of g set back to the time of the coefficient
    CHECK(g.get_time() == 1.);
  }
}
