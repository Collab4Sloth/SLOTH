/**
 * @file test_VariationalCoefficients.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests of the generated variational coefficients f(u, grad u)
 * @version 0.1
 * @date 2026-11-21
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
//
// The structures <class_name>VariationalEval are checked:
//   - against exact derivatives for simple energies,
//   - against centred finite differences of f (first derivatives) and of its first derivatives
//     (second derivatives) for every coefficient,
// and the classes through FunctionCoefficient (eval_variational_batch, point by point
// evaluation not available).

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "Coefficients/VariationalEval.hpp"
#include "sloth_unit.hpp"

// Generated at configure time from Coefficients/CoefficientExample.json:
//   EnergyIsotropic     = kappa/2 |grad phi|^2 + A phi^2 (1 - phi)^2
//   EnergyKappa         = (k0 + k1 phi)/2 |grad phi|^2
//   EnergySemiImplicit  = (1 + c^2)/2 grad u . grad ue + k/2 u |grad c|^2 + 1e-3 T X u,
//                              ue explicit, c auxiliary, (X, Y) coordinates, T time
//   EnergyAnisotropic2D = 1/2 (W0 (1 + e4 cos 4 theta))^2 |grad phi|^2 (PFHub benchmark 3)
//   EnergyAnisotropic3D = 1/2 (W0 (1 - 3 e4 + 4 e4 (nx^4 + ny^4 + nz^4)))^2 |grad phi|^2
//   EnergyTwoVariables  = 1/2 (k11 |grad a|^2 + 2 k12 grad a . grad b + k22 |grad b|^2)
//                              + a^2 b,  a in [0, 1]
#include "VariationalCases.hpp"

using Catch::Matchers::WithinAbs;

namespace {
constexpr double tol = 1.e-12;
// Finite differences: relative step and tolerance
constexpr double h = 1.e-5;
constexpr double tol_fd = 1.e-6;

/// @brief Inputs of a structure at one point
struct Point {
  std::vector<double> values, explicit_values, auxiliary_values;
  std::vector<double> gradients, explicit_gradients, auxiliary_gradients;
  std::vector<double> coordinates;
  double time = 0.;
};

/// @brief Outputs of a structure at one point
struct Values {
  double F = 0.;
  std::vector<double> df_dvar, df_dgrad, d2f_dvar2, d2f_dvar_dgrad, d2f_dgrad2;
};

/// @brief All the outputs of the structure eval at the point p (host)
template <class EVAL>
Values evaluate(const EVAL& eval, const Point& p) {
  constexpr int n = EVAL::nb_variables;
  constexpr int nd = n * EVAL::dimension;
  Values v;
  v.df_dvar.assign(n, 0.);
  v.df_dgrad.assign(nd, 0.);
  v.d2f_dvar2.assign(n * n, 0.);
  v.d2f_dvar_dgrad.assign(n * nd, 0.);
  v.d2f_dgrad2.assign(nd * nd, 0.);
  VariationalOutputs out;
  out.F = &v.F;
  out.df_dvar = v.df_dvar.data();
  out.df_dgrad = v.df_dgrad.data();
  out.d2f_dvar2 = v.d2f_dvar2.data();
  out.d2f_dvar_dgrad = v.d2f_dvar_dgrad.data();
  out.d2f_dgrad2 = v.d2f_dgrad2.data();
  eval(p.values.data(), p.explicit_values.data(), p.auxiliary_values.data(), p.gradients.data(),
       p.explicit_gradients.data(), p.auxiliary_gradients.data(), p.coordinates.data(), p.time,
       out);
  return v;
}

void check_close(const double value, const double expected, const double tolerance = tol) {
  CHECK_THAT(value, WithinAbs(expected, tolerance * (1. + std::abs(expected))));
}

/**
 * @brief Checks the derivatives given by the structure eval at the point p against centred
 * finite differences: of f for the first derivatives, of the first derivatives for the second
 * ones. Checks also the symmetry of the second derivatives.
 */
template <class EVAL>
void check_finite_differences(const EVAL& eval, const Point& p) {
  constexpr int n = EVAL::nb_variables;
  constexpr int nd = n * EVAL::dimension;
  const Values v = evaluate(eval, p);

  // Step of the components of the gradients: relative to the norm of the gradients (the
  // derivatives of a normal n = grad u / |grad u| vary on the scale |grad u|)
  double norm = 0.;
  for (const double g : p.gradients) norm += g * g;
  const double gradient_step = h * std::max(std::sqrt(norm), 1.e-3);

  // Derivatives with respect to the input x (a value or a component of a gradient)
  const auto difference = [&](std::vector<double> Point::*input, const int index,
                              const double step) {
    Point plus = p, minus = p;
    (plus.*input)[index] += step;
    (minus.*input)[index] -= step;
    const Values vp = evaluate(eval, plus), vm = evaluate(eval, minus);
    Values d;
    d.F = (vp.F - vm.F) / (2. * step);
    for (int j = 0; j < n; ++j) {
      d.df_dvar.push_back((vp.df_dvar[j] - vm.df_dvar[j]) / (2. * step));
    }
    for (int r = 0; r < nd; ++r) {
      d.df_dgrad.push_back((vp.df_dgrad[r] - vm.df_dgrad[r]) / (2. * step));
    }
    return d;
  };

  for (int i = 0; i < n; ++i) {
    CAPTURE(i);
    const Values d = difference(&Point::values, i, h * std::max(std::abs(p.values[i]), 1.));
    check_close(v.df_dvar[i], d.F, tol_fd);
    for (int j = 0; j < n; ++j) {
      CAPTURE(j);
      check_close(v.d2f_dvar2[i * n + j], d.df_dvar[j], tol_fd);
      check_close(v.d2f_dvar2[i * n + j], v.d2f_dvar2[j * n + i]);
    }
    for (int r = 0; r < nd; ++r) {
      CAPTURE(r);
      check_close(v.d2f_dvar_dgrad[i * nd + r], d.df_dgrad[r], tol_fd);
    }
  }
  for (int r = 0; r < nd; ++r) {
    CAPTURE(r);
    const Values d = difference(&Point::gradients, r, gradient_step);
    check_close(v.df_dgrad[r], d.F, tol_fd);
    for (int i = 0; i < n; ++i) {
      CAPTURE(i);
      check_close(v.d2f_dvar_dgrad[i * nd + r], d.df_dvar[i], tol_fd);
    }
    for (int s = 0; s < nd; ++s) {
      CAPTURE(s);
      check_close(v.d2f_dgrad2[r * nd + s], d.df_dgrad[s], tol_fd);
      check_close(v.d2f_dgrad2[r * nd + s], v.d2f_dgrad2[s * nd + r]);
    }
  }
}

/// @brief A point of phi with gradient (gx, gy)
Point phi_point(const double phi, const double gx, const double gy) {
  Point p;
  p.values = {phi};
  p.gradients = {gx, gy};
  return p;
}
}  // namespace

TEST_CASE("Variational coefficients: isotropic energy (exact derivatives)",
          "[Coefficients][variational]") {
  constexpr double kappa = 2.0, A = 1.5, phi = 0.3, gx = 0.7, gy = -1.2;
  const EnergyIsotropicVariationalEval eval;
  STATIC_REQUIRE(EnergyIsotropicVariationalEval::nb_variables == 1);
  STATIC_REQUIRE(EnergyIsotropicVariationalEval::dimension == 2);
  const Values v = evaluate(eval, phi_point(phi, gx, gy));

  const double W = phi * phi * (1. - phi) * (1. - phi);
  check_close(v.F, 0.5 * kappa * (gx * gx + gy * gy) + A * W);
  check_close(v.df_dvar[0], A * (2. * phi * (1. - phi) * (1. - phi) - 2. * phi * phi * (1. - phi)));
  check_close(v.df_dgrad[0], kappa * gx);
  check_close(v.df_dgrad[1], kappa * gy);
  check_close(v.d2f_dvar2[0], A * (2. - 12. * phi + 12. * phi * phi));
  check_close(v.d2f_dvar_dgrad[0], 0.);
  check_close(v.d2f_dvar_dgrad[1], 0.);
  check_close(v.d2f_dgrad2[0], kappa);
  check_close(v.d2f_dgrad2[1], 0.);
  check_close(v.d2f_dgrad2[2], 0.);
  check_close(v.d2f_dgrad2[3], kappa);
}

TEST_CASE("Variational coefficients: coefficient depending on the variable (cross terms)",
          "[Coefficients][variational]") {
  constexpr double k0 = 1.0, k1 = 0.5, phi = 0.4, gx = 0.3, gy = 0.8;
  const Values v = evaluate(EnergyKappaVariationalEval{}, phi_point(phi, gx, gy));
  check_close(v.df_dvar[0], 0.5 * k1 * (gx * gx + gy * gy));
  check_close(v.df_dgrad[0], (k0 + k1 * phi) * gx);
  check_close(v.d2f_dvar_dgrad[0], k1 * gx);
  check_close(v.d2f_dvar_dgrad[1], k1 * gy);
  check_close(v.d2f_dgrad2[0], k0 + k1 * phi);
}

TEST_CASE("Variational coefficients: explicit and auxiliary gradients, coordinates and time",
          "[Coefficients][variational]") {
  STATIC_REQUIRE(EnergySemiImplicitVariationalEval::uses_explicit_gradients);
  STATIC_REQUIRE(EnergySemiImplicitVariationalEval::uses_auxiliary_gradients);
  constexpr double k = 3.0, u = 0.6, c = 0.5, X = 2.0, T = 10.;
  Point p;
  p.values = {u};
  p.explicit_values = {0.55};
  p.auxiliary_values = {c};
  p.gradients = {0.2, -0.4};
  p.explicit_gradients = {0.25, -0.35};
  p.auxiliary_gradients = {1.5, 0.5};
  p.coordinates = {X, 7.0};
  p.time = T;
  const Values v = evaluate(EnergySemiImplicitVariationalEval{}, p);

  // Gradients of the explicit and auxiliary variables are data: not differentiated
  const double grad_c2 = 1.5 * 1.5 + 0.5 * 0.5;
  check_close(v.F, 0.5 * (1. + c * c) * (0.2 * 0.25 + 0.4 * 0.35) + 0.5 * k * u * grad_c2 +
                       1.e-3 * T * X * u);
  check_close(v.df_dvar[0], 0.5 * k * grad_c2 + 1.e-3 * T * X);
  check_close(v.df_dgrad[0], 0.5 * (1. + c * c) * 0.25);
  check_close(v.df_dgrad[1], 0.5 * (1. + c * c) * -0.35);
  for (const double value : v.d2f_dgrad2) check_close(value, 0.);
  check_finite_differences(EnergySemiImplicitVariationalEval{}, p);
}

TEST_CASE("Variational coefficients: two variables with a cross gradient term and bounds",
          "[Coefficients][variational]") {
  constexpr double k11 = 2.0, k12 = 0.5, k22 = 3.0;
  Point p;
  p.values = {0.3, 0.8};
  p.gradients = {0.1, 0.2, -0.3, 0.4};  // grad a, grad b
  const Values v = evaluate(EnergyTwoVariablesVariationalEval{}, p);

  // d2f/dp2 = K x I_2, K = [[k11, k12], [k12, k22]]
  const std::vector<double> K{k11, k12, k12, k22};
  for (int r = 0; r < 4; ++r) {
    for (int s = 0; s < 4; ++s) {
      CAPTURE(r, s);
      const double expected = (r % 2 == s % 2) ? K[(r / 2) * 2 + s / 2] : 0.;
      check_close(v.d2f_dgrad2[r * 4 + s], expected);
    }
  }
  check_finite_differences(EnergyTwoVariablesVariationalEval{}, p);

  // a clamped to [0, 1]: f and its derivatives evaluated at a = 1
  Point outside = p, at_bound = p;
  outside.values[0] = 1.5;
  at_bound.values[0] = 1.0;
  check_close(evaluate(EnergyTwoVariablesVariationalEval{}, outside).F,
              evaluate(EnergyTwoVariablesVariationalEval{}, at_bound).F);
}

TEST_CASE("Variational coefficients: anisotropic energies (finite differences)",
          "[Coefficients][variational]") {
  SECTION("four-fold anisotropy in 2D") {
    for (const auto& p :
         {phi_point(0.5, 0.7, 0.2), phi_point(0.5, -0.3, 1.1), phi_point(0.5, 1.e-3, -2.e-3)}) {
      check_finite_differences(EnergyAnisotropic2DVariationalEval{}, p);
    }
  }
  SECTION("cubic anisotropy in 3D") {
    Point p;
    p.values = {0.5};
    p.gradients = {0.3, -0.6, 0.9};
    check_finite_differences(EnergyAnisotropic3DVariationalEval{}, p);
  }
  SECTION("finite values at zero gradient (regularized normal)") {
    const Values v = evaluate(EnergyAnisotropic2DVariationalEval{}, phi_point(0.5, 0., 0.));
    CHECK(std::isfinite(v.F));
    for (const double value : v.df_dgrad) CHECK(std::isfinite(value));
    for (const double value : v.d2f_dgrad2) CHECK(std::isfinite(value));
  }
}

TEST_CASE("Variational coefficients: class and evaluation by batch",
          "[Coefficients][variational]") {
  EnergySemiImplicit coef;
  FunctionCoefficient& base = coef;
  CHECK(base.is_variational());
  CHECK(base.variational_dimension() == 2);
  CHECK(base.uses_explicit_gradients());
  CHECK(base.uses_auxiliary_gradients());
  CHECK(base.uses_coordinates());
  CHECK(base.batch_nb_variables() == 0);

  SECTION("point by point evaluation not available") {
    const std::vector<double> x{0.5};
    CHECK_THROWS_AS(base.eval_f(x), std::runtime_error);
    CHECK_THROWS_AS(base.eval_gradient(0, x), std::runtime_error);
    CHECK_THROWS_AS(base.eval_hessian(0, 0, x), std::runtime_error);
  }

  SECTION("same values as the structure at each point") {
    constexpr int npts = 5;
    constexpr double time = 4.;
    base.set_time(time);
    mfem::Vector values(npts), explicit_values(npts), auxiliary_values(npts);
    mfem::Vector gradients(2 * npts), explicit_gradients(2 * npts), auxiliary_gradients(2 * npts);
    mfem::Vector coordinates(2 * npts);
    for (int q = 0; q < npts; ++q) {
      values(q) = 0.1 * q;
      explicit_values(q) = 0.2 + 0.1 * q;
      auxiliary_values(q) = 1. - 0.1 * q;
      for (int k = 0; k < 2; ++k) {
        gradients(2 * q + k) = 0.3 * q - 0.2 * k;
        explicit_gradients(2 * q + k) = 0.1 + 0.05 * q * (k + 1);
        auxiliary_gradients(2 * q + k) = -0.4 + 0.2 * q + k;
        coordinates(2 * q + k) = 1. + q + 2. * k;
      }
    }
    // Only F and df/dgrad requested: the other outputs are not computed
    mfem::Vector F(npts), df_dgrad(2 * npts);
    VariationalBatch outputs;
    outputs.F = &F;
    outputs.df_dgrad = &df_dgrad;
    REQUIRE(base.eval_variational_batch(npts, values, explicit_values, auxiliary_values, gradients,
                                        explicit_gradients, auxiliary_gradients, coordinates,
                                        outputs));
    F.HostRead();
    df_dgrad.HostRead();
    for (int q = 0; q < npts; ++q) {
      CAPTURE(q);
      Point p;
      p.values = {values(q)};
      p.explicit_values = {explicit_values(q)};
      p.auxiliary_values = {auxiliary_values(q)};
      p.gradients = {gradients(2 * q), gradients(2 * q + 1)};
      p.explicit_gradients = {explicit_gradients(2 * q), explicit_gradients(2 * q + 1)};
      p.auxiliary_gradients = {auxiliary_gradients(2 * q), auxiliary_gradients(2 * q + 1)};
      p.coordinates = {coordinates(2 * q), coordinates(2 * q + 1)};
      p.time = time;
      const Values v = evaluate(coef.get_eval(), p);
      check_close(F(q), v.F);
      check_close(df_dgrad(2 * q), v.df_dgrad[0]);
      check_close(df_dgrad(2 * q + 1), v.df_dgrad[1]);
    }
  }
}
