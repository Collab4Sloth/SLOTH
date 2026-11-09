/**
 * @file test_GeneratedCoefficients.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests for coefficients generated from JSON by GenerateCoefficient.py
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
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <span>
#include <vector>

#include "Coefficients/Coefficients.hpp"
#include "Coefficients/EvalAtPoints.hpp"
#include "Coefficients/ProductFunctionCoefficients.hpp"
#include "Coefficients/SumFunctionCoefficients.hpp"
#include "sloth_unit.hpp"

// Generated at configure time from Coefficients/CoefficientExample.json:
//   FunctionA = a x + b y + c z,  (a, b, c) = (2, 30, 0.1)
//   FunctionB = 2 x^2 + 30 x y + 0.1 x y z
//   FunctionC to FunctionH: see the evaluation structures and the bounds below
#include "EvalCases.hpp"
#include "Function.hpp"
#include "OtherFunction.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
constexpr double tol = 1.e-10;
// Evaluation point (x, y, z)
constexpr std::array<double, 3> X{2., 3., 1.};

// Exact values at X, computed by hand
constexpr double A = 94.1;   // 2*2 + 30*3 + 0.1*1
constexpr double B = 188.6;  // 2*4 + 30*6 + 0.1*6
const std::vector<double> gradA{2., 30., 0.1};
const std::vector<double> gradB{98.3, 60.2, 0.6};
const std::vector<double> hessB{4.0,  30.1, 0.3,  //
                                30.1, 0.0,  0.2,  //
                                0.3,  0.2,  0.0};

/// @brief Checks the value, gradient and Hessian of c at X
void check_at_X(Coefficient& c, double value, const std::vector<double>& grad,
                const std::vector<double>& hess) {
  CHECK_THAT(c.compute(X), WithinRel(value, tol));
  for (int i = 0; i < 3; ++i) {
    CAPTURE(i);
    CHECK_THAT(c.compute_gradient(i, X), WithinRel(grad[i], tol));
    for (int j = 0; j < 3; ++j) {
      CAPTURE(j);
      const double expected = hess[i * 3 + j];
      CHECK_THAT(c.compute_hessian(i, j, X), WithinAbs(expected, tol * (1. + std::abs(expected))));
    }
  }
}
}  // namespace

TEST_CASE("Generated coefficients: linear function (FunctionA)", "[Coefficients][generated]") {
  Coefficient c(Glossary::Temperature, Scheme::Implicit, FunctionA());
  check_at_X(c, A, gradA, std::vector<double>(9, 0.));
}

TEST_CASE("Generated coefficients: quadratic function (FunctionB)", "[Coefficients][generated]") {
  Coefficient c(Glossary::Temperature, Scheme::Implicit, FunctionB());
  check_at_X(c, B, gradB, hessB);
}

TEST_CASE("Generated coefficients: in a Coefficients collection", "[Coefficients][generated]") {
  const Coefficient ca(Glossary::Temperature, Scheme::Implicit, FunctionA());
  const Coefficient cb(Glossary::Temperature, Scheme::Implicit, FunctionB());
  Coefficients coefs(ca, cb);
  coefs.add(cb);
  REQUIRE(coefs.size() == 3);
  check_at_X(coefs[1], B, gradB, hessB);
}

TEST_CASE("Generated coefficients: sum B + B", "[Coefficients][generated]") {
  FunctionB fb;
  Coefficient c(Glossary::Temperature, Scheme::Implicit, SumCoefficient({&fb, &fb}));
  std::vector<double> grad2, hess2;
  for (double g : gradB) grad2.push_back(2. * g);
  for (double h : hessB) hess2.push_back(2. * h);
  check_at_X(c, 2. * B, grad2, hess2);
}

TEST_CASE("Generated coefficients: product A * B", "[Coefficients][generated]") {
  FunctionA fa;
  FunctionB fb;
  Coefficient c(Glossary::Temperature, Scheme::Implicit, ProductCoefficient({&fa, &fb}));

  // (A B)' = A' B + A B',  (A B)'' = A'' B + A' B'^T + B' A'^T + A B''   (A'' = 0)
  std::vector<double> grad(3), hess(9);
  for (int i = 0; i < 3; ++i) {
    grad[i] = gradA[i] * B + A * gradB[i];
    for (int j = 0; j < 3; ++j) {
      hess[i * 3 + j] = gradA[i] * gradB[j] + gradB[i] * gradA[j] + A * hessB[i * 3 + j];
    }
  }
  check_at_X(c, A * B, grad, hess);

  SECTION("sum of products: P + P with P = A * B") {
    auto p = ProductCoefficient({&fa, &fb});
    Coefficient s(Glossary::Temperature, Scheme::Implicit, SumCoefficient({&p, &p}));
    std::vector<double> grad2, hess2;
    for (double g : grad) grad2.push_back(2. * g);
    for (double h : hess) hess2.push_back(2. * h);
    check_at_X(s, 2. * A * B, grad2, hess2);
  }
}

// Sums and products of generated coefficients evaluated by batch (eval_batch): same values,
// gradients and Hessians as the point by point evaluation (eval_f, eval_gradient, eval_hessian)
namespace {
/// @brief Compares eval_batch of c with its point by point evaluation at npts points (x, y, z)
void check_batch(FunctionCoefficient& c, const bool with_F, const bool with_gradient,
                 const bool with_hessian) {
  constexpr int npts = 4;
  constexpr int n = 3;
  REQUIRE(c.batch_nb_variables() == n);
  mfem::Vector values(n * npts), none;
  for (int q = 0; q < npts; ++q) {
    values(n * q) = 0.5 + q;
    values(n * q + 1) = -1. + 0.7 * q;
    values(n * q + 2) = 2. - 0.3 * q;
  }
  mfem::Vector F(npts), gradient(n * npts), hessian(n * n * npts);
  REQUIRE(c.eval_batch(npts, values, none, none, none, with_F ? &F : nullptr,
                       with_gradient ? &gradient : nullptr, with_hessian ? &hessian : nullptr));
  F.HostRead();
  gradient.HostRead();
  hessian.HostRead();
  values.HostRead();
  for (int q = 0; q < npts; ++q) {
    CAPTURE(q);
    const std::span<const double> x(values.GetData() + n * q, n);
    if (with_F) {
      const double expected = c.eval_f(x);
      CHECK_THAT(F(q), WithinAbs(expected, tol * (1. + std::abs(expected))));
    }
    for (int i = 0; i < n; ++i) {
      CAPTURE(i);
      if (with_gradient) {
        const double expected = c.eval_gradient(i, x);
        CHECK_THAT(gradient(n * q + i), WithinAbs(expected, tol * (1. + std::abs(expected))));
      }
      for (int j = 0; j < n && with_hessian; ++j) {
        CAPTURE(j);
        const double expected = c.eval_hessian(i, j, x);
        CHECK_THAT(hessian((q * n + i) * n + j),
                   WithinAbs(expected, tol * (1. + std::abs(expected))));
      }
    }
  }
}
}  // namespace

TEST_CASE("Sum and product of coefficients: evaluation by batch",
          "[Coefficients][generated][eval]") {
  FunctionA fa;
  FunctionB fb;
  SECTION("sum A + B") {
    SumCoefficient c({&fa, &fb});
    check_batch(c, true, true, true);
  }
  SECTION("product A * B") {
    ProductCoefficient c({&fa, &fb});
    check_batch(c, true, true, true);
  }
  SECTION("product A * B * B") {
    ProductCoefficient c({&fa, &fb, &fb});
    check_batch(c, true, true, true);
  }
  SECTION("product A * B: gradient only, Hessian only") {
    ProductCoefficient c({&fa, &fb, &fb});
    check_batch(c, false, true, false);
    check_batch(c, false, false, true);
  }
  SECTION("nested: (A * B) + (A * B) and (A + B) * B") {
    ProductCoefficient p({&fa, &fb});
    SumCoefficient s({&p, &p});
    check_batch(s, true, true, true);
    SumCoefficient a_plus_b({&fa, &fb});
    ProductCoefficient q({&a_plus_b, &fb});
    check_batch(q, true, true, true);
  }
  SECTION("terms with different numbers of variables: point by point evaluation") {
    FunctionC fc;  // one variable
    ProductCoefficient c({&fa, &fc});
    CHECK(c.batch_nb_variables() == 0);
    mfem::Vector values(3), none, F(1);
    CHECK_FALSE(c.eval_batch(1, values, none, none, none, &F, nullptr, nullptr));
  }
}

//////////////////////////////////////////////////////////////////////////////
// Evaluation structures <class_name>Eval (host and device)
//////////////////////////////////////////////////////////////////////////////
//
// Generated at configure time from Coefficients/CoefficientExample.json (EvalCases.hpp):
//   FunctionC = a c^2 (1 - c) + R T c (eta1^2 + eta2^2),  a = 4, R = Physical::R, T = time
//   FunctionD = x^2 y + 3 x xe ye,  (xe, ye) explicit variables
//   FunctionE = v (1 - ((X - 1)/2)^2) u^2,  v = 0.5, (X, Y) coordinates
//   FunctionF = x^2 on [0, 1], 2 x - 1 on ]1, 2], 3 elsewhere

namespace {
/// @brief Values of a coefficient at one point
struct EvalOutput {
  double F = 0.;
  std::vector<double> gradient;
  std::vector<double> hessian;
};

/// @brief Evaluates the structure eval at one point
template <class EVAL>
EvalOutput eval_at_point(const EVAL& eval, const std::vector<double>& values,
                         const std::vector<double>& explicit_values,
                         const std::vector<double>& auxiliary_values,
                         const std::vector<double>& coordinates, const double time) {
  constexpr int n = EVAL::nb_variables;
  EvalOutput out;
  out.gradient.assign(n, 0.);
  out.hessian.assign(n * n, 0.);
  eval(values.data(), explicit_values.data(), auxiliary_values.data(), coordinates.data(), time,
       &out.F, out.gradient.data(), out.hessian.data());
  return out;
}

/// @brief Checks that the structure eval gives the same F, gradient and Hessian as the class coef
template <class EVAL>
void check_eval_as_class(FunctionCoefficient& coef, const EVAL& eval,
                         const std::vector<double>& values,
                         const std::vector<double>& explicit_values,
                         const std::vector<double>& auxiliary_values,
                         const std::vector<double>& coordinates, const double time) {
  constexpr int n = EVAL::nb_variables;
  REQUIRE(static_cast<int>(values.size()) == n);
  coef.set_time(time);
  if (!coordinates.empty()) {
    coef.set_coordinates(coordinates);
  }
  const std::span<const double> u(values);
  const std::span<const double> u_exp(explicit_values);
  const std::span<const double> aux(auxiliary_values);

  const EvalOutput out =
      eval_at_point(eval, values, explicit_values, auxiliary_values, coordinates, time);

  const double f = coef.eval_f(u, u_exp, aux);
  CHECK_THAT(out.F, WithinAbs(f, tol * (1. + std::abs(f))));
  for (int i = 0; i < n; ++i) {
    CAPTURE(i);
    const double g = coef.eval_gradient(i, u, u_exp, aux);
    CHECK_THAT(out.gradient[i], WithinAbs(g, tol * (1. + std::abs(g))));
    for (int j = 0; j < n; ++j) {
      CAPTURE(j);
      const double h = coef.eval_hessian(i, j, u, u_exp, aux);
      CHECK_THAT(out.hessian[i * n + j], WithinAbs(h, tol * (1. + std::abs(h))));
    }
  }
}
}  // namespace

TEST_CASE("Evaluation structures: exact values (FunctionA, FunctionB)",
          "[Coefficients][generated][eval]") {
  const std::vector<double> x(X.begin(), X.end());

  const EvalOutput a = eval_at_point(FunctionAEval{}, x, {}, {}, {}, 0.);
  CHECK_THAT(a.F, WithinRel(A, tol));
  const EvalOutput b = eval_at_point(FunctionBEval{}, x, {}, {}, {}, 0.);
  CHECK_THAT(b.F, WithinRel(B, tol));
  for (int i = 0; i < 3; ++i) {
    CAPTURE(i);
    CHECK_THAT(a.gradient[i], WithinRel(gradA[i], tol));
    CHECK_THAT(b.gradient[i], WithinRel(gradB[i], tol));
    for (int j = 0; j < 3; ++j) {
      CAPTURE(j);
      const double expected = hessB[i * 3 + j];
      CHECK_THAT(a.hessian[i * 3 + j], WithinAbs(0., tol));
      CHECK_THAT(b.hessian[i * 3 + j], WithinAbs(expected, tol * (1. + std::abs(expected))));
    }
  }
}

TEST_CASE("Evaluation structures: prefactor", "[Coefficients][generated][eval]") {
  const std::vector<double> x(X.begin(), X.end());
  FunctionB coef(3.);
  FunctionBEval eval;
  eval.prefactor_ = 3.;
  check_eval_as_class(coef, eval, x, {}, {}, {}, 0.);
}

TEST_CASE("Evaluation structures: same values as the classes", "[Coefficients][generated][eval]") {
  SECTION("auxiliary variables, constants and time (FunctionC)") {
    FunctionC coef;
    for (const double c : {0.1, 0.5, 0.9}) {
      CAPTURE(c);
      check_eval_as_class(coef, FunctionCEval{}, {c}, {}, {0.3, 0.7}, {}, 1200.);
    }
  }
  SECTION("explicit variables (FunctionD)") {
    FunctionD coef;
    check_eval_as_class(coef, FunctionDEval{}, {2., 3.}, {0.5, -1.5}, {}, {}, 0.);
    check_eval_as_class(coef, FunctionDEval{}, {-1., 0.25}, {4., 2.}, {}, {}, 0.);
  }
  SECTION("coordinates (FunctionE)") {
    FunctionE coef;
    check_eval_as_class(coef, FunctionEEval{}, {0.8}, {}, {}, {0.4, 2.}, 0.);
    check_eval_as_class(coef, FunctionEEval{}, {-2.}, {}, {}, {1.7, -1.}, 0.);
  }
  SECTION("conditional expression, on each piece and at the bounds (FunctionF)") {
    FunctionF coef;
    for (const double x : {-1., 0., 0.5, 1., 1.5, 2., 2.5}) {
      CAPTURE(x);
      check_eval_as_class(coef, FunctionFEval{}, {x}, {}, {}, {}, 0.);
    }
  }
}

TEST_CASE("Evaluation structures: evaluation at a set of points (eval_at_points)",
          "[Coefficients][generated][eval]") {
  // FunctionC at npts points: c_q = 0.1 + 0.2 q, eta_q = (0.1 q, 1 - 0.1 q)
  constexpr int npts = 5;
  constexpr double time = 300.;
  mfem::Vector values(npts), auxiliary(2 * npts), none;
  for (int q = 0; q < npts; ++q) {
    values(q) = 0.1 + 0.2 * q;
    auxiliary(2 * q) = 0.1 * q;
    auxiliary(2 * q + 1) = 1. - 0.1 * q;
  }
  const FunctionCEval eval;

  SECTION("F, gradient and Hessian") {
    mfem::Vector F(npts), gradient(npts), hessian(npts);
    eval_at_points(eval, npts, time, values, none, auxiliary, none, &F, &gradient, &hessian);
    F.HostRead();
    gradient.HostRead();
    hessian.HostRead();
    for (int q = 0; q < npts; ++q) {
      CAPTURE(q);
      const EvalOutput expected =
          eval_at_point(eval, {values(q)}, {}, {auxiliary(2 * q), auxiliary(2 * q + 1)}, {}, time);
      CHECK_THAT(F(q), WithinRel(expected.F, tol));
      CHECK_THAT(gradient(q), WithinRel(expected.gradient[0], tol));
      const double h = expected.hessian[0];
      CHECK_THAT(hessian(q), WithinAbs(h, tol * (1. + std::abs(h))));
    }
  }

  SECTION("gradient only") {
    mfem::Vector gradient(npts);
    eval_at_points(eval, npts, time, values, none, auxiliary, none, nullptr, &gradient);
    gradient.HostRead();
    for (int q = 0; q < npts; ++q) {
      CAPTURE(q);
      const EvalOutput expected =
          eval_at_point(eval, {values(q)}, {}, {auxiliary(2 * q), auxiliary(2 * q + 1)}, {}, time);
      CHECK_THAT(gradient(q), WithinRel(expected.gradient[0], tol));
    }
  }

  SECTION("inputs with more components per point than read by the structure") {
    // 3 values and 4 auxiliary fields per point: FunctionC reads the first 1 and the first 2
    mfem::Vector values3(3 * npts), auxiliary4(4 * npts), F(npts), gradient(npts);
    for (int q = 0; q < npts; ++q) {
      values3(3 * q) = values(q);
      values3(3 * q + 1) = 100.;
      values3(3 * q + 2) = -100.;
      auxiliary4(4 * q) = auxiliary(2 * q);
      auxiliary4(4 * q + 1) = auxiliary(2 * q + 1);
      auxiliary4(4 * q + 2) = 50.;
      auxiliary4(4 * q + 3) = -50.;
    }
    eval_at_points(eval, npts, time, values3, none, auxiliary4, none, &F, &gradient);
    F.HostRead();
    gradient.HostRead();
    values.HostRead();
    auxiliary.HostRead();
    for (int q = 0; q < npts; ++q) {
      CAPTURE(q);
      const EvalOutput expected =
          eval_at_point(eval, {values(q)}, {}, {auxiliary(2 * q), auxiliary(2 * q + 1)}, {}, time);
      CHECK_THAT(F(q), WithinRel(expected.F, tol));
      CHECK_THAT(gradient(q), WithinRel(expected.gradient[0], tol));
    }
  }

  SECTION("two variables: gradient and Hessian of each point are contiguous (FunctionD)") {
    constexpr int n = FunctionDEval::nb_variables;
    mfem::Vector u(n * npts), u_exp(n * npts), gradient(n * npts), hessian(n * n * npts);
    for (int q = 0; q < npts; ++q) {
      u(n * q) = 1. + q;
      u(n * q + 1) = 0.5 * q;
      u_exp(n * q) = -1. + 0.3 * q;
      u_exp(n * q + 1) = 2.;
    }
    eval_at_points(FunctionDEval{}, npts, 0., u, u_exp, none, none, nullptr, &gradient, &hessian);
    gradient.HostRead();
    hessian.HostRead();
    for (int q = 0; q < npts; ++q) {
      CAPTURE(q);
      const EvalOutput expected = eval_at_point(FunctionDEval{}, {u(n * q), u(n * q + 1)},
                                                {u_exp(n * q), u_exp(n * q + 1)}, {}, {}, 0.);
      for (int i = 0; i < n; ++i) {
        CAPTURE(i);
        CHECK_THAT(gradient(n * q + i),
                   WithinAbs(expected.gradient[i], tol * (1. + std::abs(expected.gradient[i]))));
        for (int j = 0; j < n; ++j) {
          CAPTURE(j);
          const double h = expected.hessian[i * n + j];
          CHECK_THAT(hessian(n * n * q + i * n + j), WithinAbs(h, tol * (1. + std::abs(h))));
        }
      }
    }
  }

  SECTION("inconsistent sizes are rejected") {
    mfem::Vector F(npts + 1);
    SLOTH_CHECK_MFEM_ERROR(eval_at_points(eval, npts, time, values, none, auxiliary, none, &F));
    mfem::Vector F_ok(npts);
    SLOTH_CHECK_MFEM_ERROR(eval_at_points(eval, npts, time, values, none, none, none, &F_ok));
  }
}

//////////////////////////////////////////////////////////////////////////////
// Bounds of variables ("bounds" attribute of the JSON file)
//////////////////////////////////////////////////////////////////////////////
//
// Generated at configure time from Coefficients/CoefficientExample.json (EvalCases.hpp):
//   FunctionG = x^2 y + 2 y eta1 eta2 + x eta2,
//               x in [lo, hi] = [-1, 2], y in [0, +inf[, eta1, eta2 in ]-inf, 1]
//   FunctionH = x^2 on [0, 1], 2 x - 1 on ]1, 2], 100 elsewhere, x in [0, 2]
// A bounded variable is clamped to its bounds: F, gradient and Hessian are evaluated at the
// clamped values.

namespace {
/// @brief Exact F, gradient and Hessian of FunctionG at (x, y, eta1, eta2) (inside the bounds)
EvalOutput exact_G(const double x, const double y, const double eta1, const double eta2) {
  EvalOutput out;
  out.F = x * x * y + 2. * y * eta1 * eta2 + x * eta2;
  out.gradient = {2. * x * y + eta2, x * x + 2. * eta1 * eta2};
  out.hessian = {2. * y, 2. * x, 2. * x, 0.};
  return out;
}

void check_output(const EvalOutput& out, const EvalOutput& expected) {
  CHECK_THAT(out.F, WithinAbs(expected.F, tol * (1. + std::abs(expected.F))));
  REQUIRE(out.gradient.size() == expected.gradient.size());
  REQUIRE(out.hessian.size() == expected.hessian.size());
  for (std::size_t i = 0; i < out.gradient.size(); ++i) {
    CAPTURE(i);
    const double g = expected.gradient[i];
    CHECK_THAT(out.gradient[i], WithinAbs(g, tol * (1. + std::abs(g))));
  }
  for (std::size_t k = 0; k < out.hessian.size(); ++k) {
    CAPTURE(k);
    const double h = expected.hessian[k];
    CHECK_THAT(out.hessian[k], WithinAbs(h, tol * (1. + std::abs(h))));
  }
}
}  // namespace

TEST_CASE("Bounds of variables: values inside the bounds are not modified",
          "[Coefficients][generated][eval][bounds]") {
  const EvalOutput out = eval_at_point(FunctionGEval{}, {0.5, 0.3}, {}, {0.2, 0.4}, {}, 0.);
  check_output(out, exact_G(0.5, 0.3, 0.2, 0.4));

  SECTION("bounds reached exactly") {
    check_output(eval_at_point(FunctionGEval{}, {2., 0.}, {}, {1., 1.}, {}, 0.),
                 exact_G(2., 0., 1., 1.));
    check_output(eval_at_point(FunctionGEval{}, {-1., 0.}, {}, {1., -3.}, {}, 0.),
                 exact_G(-1., 0., 1., -3.));
  }
}

TEST_CASE("Bounds of variables: values outside the bounds are clamped",
          "[Coefficients][generated][eval][bounds]") {
  SECTION("both bounds (x), lower bound only (y), upper bound only (eta1, eta2)") {
    // (x, y, eta1, eta2) = (3, -1, 1.5, 2) is clamped to (2, 0, 1, 1)
    const EvalOutput out = eval_at_point(FunctionGEval{}, {3., -1.}, {}, {1.5, 2.}, {}, 0.);
    check_output(out, exact_G(2., 0., 1., 1.));
    // Values by hand: F = 2, gradient = (1, 6), Hessian = ((0, 4), (4, 0))
    CHECK_THAT(out.F, WithinAbs(2., tol));
    CHECK_THAT(out.gradient[0], WithinAbs(1., tol));
    CHECK_THAT(out.gradient[1], WithinAbs(6., tol));
    CHECK_THAT(out.hessian[1], WithinAbs(4., tol));
  }
  SECTION("one-sided bounds: no clamping on the side without bound") {
    // x = -5 is clamped to lo = -1; y = 5 and eta = (-2, -3) have no bound on this side
    check_output(eval_at_point(FunctionGEval{}, {-5., 5.}, {}, {-2., -3.}, {}, 0.),
                 exact_G(-1., 5., -2., -3.));
  }
  SECTION("the class gives the same values as its evaluation structure") {
    FunctionG coef;
    check_eval_as_class(coef, FunctionGEval{}, {3., -1.}, {}, {1.5, 2.}, {}, 0.);
    check_eval_as_class(coef, FunctionGEval{}, {-5., 5.}, {}, {-2., -3.}, {}, 0.);
  }
  SECTION("constants used as bounds follow the members of the structure") {
    FunctionGEval eval;
    eval.hi = 4.;
    check_output(eval_at_point(eval, {3., 1.}, {}, {0.5, 0.5}, {}, 0.), exact_G(3., 1., 0.5, 0.5));
  }
}

TEST_CASE("Bounds of variables: conditional expression on a bounded variable",
          "[Coefficients][generated][eval][bounds]") {
  // The piece is selected with the clamped value: x = 3 -> 2 (second piece), x = -1 -> 0 (first
  // piece); the default expression (100) is never used
  EvalOutput expected;
  expected.F = 3.;
  expected.gradient = {2.};
  expected.hessian = {0.};
  check_output(eval_at_point(FunctionHEval{}, {3.}, {}, {}, {}, 0.), expected);

  expected.F = 0.;
  expected.gradient = {0.};
  expected.hessian = {2.};
  check_output(eval_at_point(FunctionHEval{}, {-1.}, {}, {}, {}, 0.), expected);
}

TEST_CASE("Bounds of variables: evaluation at a set of points (eval_at_points)",
          "[Coefficients][generated][eval][bounds]") {
  // Points inside and outside the bounds
  constexpr int npts = 5;
  constexpr int n = FunctionGEval::nb_variables;
  const double xs[npts] = {-3., -1., 0.5, 2., 4.};
  const double ys[npts] = {-2., 0., 0.3, 1., 3.};
  mfem::Vector u(n * npts), aux(2 * npts), F(npts), gradient(n * npts), none;
  for (int q = 0; q < npts; ++q) {
    u(n * q) = xs[q];
    u(n * q + 1) = ys[q];
    aux(2 * q) = 0.5 * q;
    aux(2 * q + 1) = 1.5 - 0.5 * q;
  }
  eval_at_points(FunctionGEval{}, npts, 0., u, none, aux, none, &F, &gradient);
  F.HostRead();
  gradient.HostRead();
  for (int q = 0; q < npts; ++q) {
    CAPTURE(q);
    const double x = std::clamp(xs[q], -1., 2.);
    const double y = std::max(ys[q], 0.);
    const double eta1 = std::min(aux(2 * q), 1.);
    const double eta2 = std::min(aux(2 * q + 1), 1.);
    const EvalOutput expected = exact_G(x, y, eta1, eta2);
    CHECK_THAT(F(q), WithinAbs(expected.F, tol * (1. + std::abs(expected.F))));
    for (int i = 0; i < n; ++i) {
      CAPTURE(i);
      const double g = expected.gradient[static_cast<std::size_t>(i)];
      CHECK_THAT(gradient(n * q + i), WithinAbs(g, tol * (1. + std::abs(g))));
    }
  }
}
