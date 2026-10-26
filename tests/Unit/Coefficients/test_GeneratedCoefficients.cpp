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
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <vector>

#include "Coefficients/Coefficients.hpp"
#include "Coefficients/ProductFunctionCoefficients.hpp"
#include "Coefficients/SumFunctionCoefficients.hpp"
#include "sloth_unit.hpp"

// Generated at configure time from Coefficients/CoefficientExample.json:
//   FunctionA = a x + b y + c z,  (a, b, c) = (2, 30, 0.1)
//   FunctionB = 2 x^2 + 30 x y + 0.1 x y z
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
