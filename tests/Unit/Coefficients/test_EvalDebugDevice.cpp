/**
 * @file test_EvalDebugDevice.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Unit tests of eval_at_points with the MFEM "debug" device
 * @version 0.1
 * @date 2026-11-07
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
// The "debug" device of MFEM runs the mfem::forall loops on the host, but manages the memory as
// with a GPU: each mfem::Vector has a host copy and a device copy, and the copy which is not up
// to date is protected. Reading or writing data on the wrong side stops the program at once.
// These tests check that eval_at_points reads its inputs and writes its outputs on the device,
// and that the results are available on the host after HostRead().
//
// The device is configured once for the whole program: these tests have their own executable.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <vector>

#include "Coefficients/EvalAtPoints.hpp"
#include "Coefficients/ProductFunctionCoefficients.hpp"
#include "Coefficients/SumFunctionCoefficients.hpp"
#include "Coefficients/VariationalEval.hpp"
#include "sloth_unit.hpp"

// Generated at configure time from Coefficients/CoefficientExample.json:
//   FunctionC = a c^2 (1 - c) + R T c (eta1^2 + eta2^2),  a = 4, R = Physical::R, T = time
//   FunctionD = x^2 y + 3 x xe ye,  (xe, ye) explicit variables
//   FunctionE = v (1 - ((X - 1)/2)^2) u^2,  v = 0.5, (X, Y) coordinates
#include "EvalCases.hpp"
// EnergySemiImplicit = (1 + c^2)/2 grad u . grad ue + k/2 u |grad c|^2 + 1e-3 T X u,  k = 3,
// ue explicit, c auxiliary, (X, Y) coordinates, T time
#include "VariationalCases.hpp"

using Catch::Matchers::WithinAbs;

namespace {
constexpr double tol = 1.e-10;
constexpr int npts = 6;

/// @brief Configures the MFEM "debug" device (once for the whole program)
void use_debug_device() {
  static mfem::Device device("debug");
  REQUIRE(mfem::Device::IsEnabled());
  REQUIRE(mfem::Device::Allows(mfem::Backend::DEBUG_DEVICE));
}

/// @brief Values of the structure eval at one point (host)
template <class EVAL>
void eval_on_host(const EVAL& eval, const std::vector<double>& values,
                  const std::vector<double>& explicit_values,
                  const std::vector<double>& auxiliary_values,
                  const std::vector<double>& coordinates, const double time, double& F,
                  std::vector<double>& gradient, std::vector<double>& hessian) {
  constexpr int n = EVAL::nb_variables;
  gradient.assign(n, 0.);
  hessian.assign(n * n, 0.);
  eval(values.data(), explicit_values.data(), auxiliary_values.data(), coordinates.data(), time, &F,
       gradient.data(), hessian.data());
}

void check_close(const double value, const double expected) {
  CHECK_THAT(value, WithinAbs(expected, tol * (1. + std::abs(expected))));
}
}  // namespace

TEST_CASE("eval_at_points on the debug device: auxiliary variables and time (FunctionC)",
          "[Coefficients][eval][device]") {
  use_debug_device();
  constexpr double time = 300.;
  mfem::Vector values(npts), auxiliary(2 * npts), none;
  for (int q = 0; q < npts; ++q) {
    values(q) = 0.05 + 0.15 * q;
    auxiliary(2 * q) = 0.1 * q;
    auxiliary(2 * q + 1) = 1. - 0.1 * q;
  }
  mfem::Vector F(npts), gradient(npts), hessian(npts);
  const FunctionCEval eval;
  eval_at_points(eval, npts, time, values, none, auxiliary, none, &F, &gradient, &hessian);

  // Results written on the device: copied back to the host before being read
  F.HostRead();
  gradient.HostRead();
  hessian.HostRead();
  values.HostRead();
  auxiliary.HostRead();
  for (int q = 0; q < npts; ++q) {
    CAPTURE(q);
    double f = 0.;
    std::vector<double> g, h;
    eval_on_host(eval, {values(q)}, {}, {auxiliary(2 * q), auxiliary(2 * q + 1)}, {}, time, f, g,
                 h);
    check_close(F(q), f);
    check_close(gradient(q), g[0]);
    check_close(hessian(q), h[0]);
  }
}

TEST_CASE("eval_at_points on the debug device: two variables and explicit variables (FunctionD)",
          "[Coefficients][eval][device]") {
  use_debug_device();
  constexpr int n = FunctionDEval::nb_variables;
  mfem::Vector u(n * npts), u_exp(n * npts), none;
  for (int q = 0; q < npts; ++q) {
    u(n * q) = 1. + q;
    u(n * q + 1) = 0.5 * q;
    u_exp(n * q) = -1. + 0.3 * q;
    u_exp(n * q + 1) = 2.;
  }

  SECTION("gradient and Hessian only") {
    mfem::Vector gradient(n * npts), hessian(n * n * npts);
    eval_at_points(FunctionDEval{}, npts, 0., u, u_exp, none, none, nullptr, &gradient, &hessian);
    gradient.HostRead();
    hessian.HostRead();
    u.HostRead();
    u_exp.HostRead();
    for (int q = 0; q < npts; ++q) {
      CAPTURE(q);
      double f = 0.;
      std::vector<double> g, h;
      eval_on_host(FunctionDEval{}, {u(n * q), u(n * q + 1)}, {u_exp(n * q), u_exp(n * q + 1)}, {},
                   {}, 0., f, g, h);
      for (int i = 0; i < n; ++i) {
        CAPTURE(i);
        check_close(gradient(n * q + i), g[static_cast<std::size_t>(i)]);
        for (int j = 0; j < n; ++j) {
          CAPTURE(j);
          check_close(hessian(n * n * q + i * n + j), h[static_cast<std::size_t>(i * n + j)]);
        }
      }
    }
  }

  SECTION("inputs modified on the device between two evaluations") {
    // u is scaled by 2 on the device: the second evaluation must read the device copy
    mfem::Vector F1(npts), F2(npts);
    eval_at_points(FunctionDEval{}, npts, 0., u, u_exp, none, none, &F1);
    u.UseDevice(true);
    double* d_u = u.ReadWrite();
    mfem::forall(n * npts, [=] MFEM_HOST_DEVICE(int k) { d_u[k] *= 2.; });
    eval_at_points(FunctionDEval{}, npts, 0., u, u_exp, none, none, &F2);
    F1.HostRead();
    F2.HostRead();
    u.HostRead();
    u_exp.HostRead();
    for (int q = 0; q < npts; ++q) {
      CAPTURE(q);
      // F = x^2 y + 3 x xe ye with (x, y) the scaled values
      const double x = u(n * q), y = u(n * q + 1);
      check_close(F2(q), x * x * y + 3. * x * u_exp(n * q) * u_exp(n * q + 1));
      // F1 with the values before scaling (x / 2, y / 2)
      check_close(F1(q), 0.125 * x * x * y + 1.5 * x * u_exp(n * q) * u_exp(n * q + 1));
    }
  }
}

TEST_CASE("eval_at_points on the debug device: coordinates (FunctionE)",
          "[Coefficients][eval][device]") {
  use_debug_device();
  constexpr int dim = FunctionEEval::nb_coordinates;
  mfem::Vector values(npts), coordinates(dim * npts), F(npts), none;
  for (int q = 0; q < npts; ++q) {
    values(q) = -1. + 0.4 * q;
    coordinates(dim * q) = 0.2 * q;
    coordinates(dim * q + 1) = 1.;
  }
  eval_at_points(FunctionEEval{}, npts, 0., values, none, none, coordinates, &F);
  F.HostRead();
  values.HostRead();
  coordinates.HostRead();
  for (int q = 0; q < npts; ++q) {
    CAPTURE(q);
    const double X = coordinates(dim * q);
    const double u = values(q);
    check_close(F(q), 0.5 * (1. - 0.25 * (X - 1.) * (X - 1.)) * u * u);
  }
}

TEST_CASE("Sum and product on the debug device: D + D and D * D (FunctionD)",
          "[Coefficients][eval][device]") {
  use_debug_device();
  constexpr int n = FunctionDEval::nb_variables;
  mfem::Vector u(n * npts), u_exp(n * npts), none;
  for (int q = 0; q < npts; ++q) {
    u(n * q) = 1. + q;
    u(n * q + 1) = 0.5 * q;
    u_exp(n * q) = -1. + 0.3 * q;
    u_exp(n * q + 1) = 2.;
  }
  FunctionD fd;
  SumCoefficient sum({&fd, &fd});
  ProductCoefficient product({&fd, &fd});
  mfem::Vector F_sum(npts), g_sum(n * npts), h_sum(n * n * npts);
  mfem::Vector F_prod(npts), g_prod(n * npts), h_prod(n * n * npts);
  REQUIRE(sum.eval_batch(npts, u, u_exp, none, none, &F_sum, &g_sum, &h_sum));
  REQUIRE(product.eval_batch(npts, u, u_exp, none, none, &F_prod, &g_prod, &h_prod));

  // Results written on the device: copied back to the host before being read
  for (mfem::Vector* v : {&F_sum, &g_sum, &h_sum, &F_prod, &g_prod, &h_prod, &u, &u_exp}) {
    v->HostRead();
  }
  for (int q = 0; q < npts; ++q) {
    CAPTURE(q);
    double f = 0.;
    std::vector<double> g, h;
    eval_on_host(FunctionDEval{}, {u(n * q), u(n * q + 1)}, {u_exp(n * q), u_exp(n * q + 1)}, {},
                 {}, 0., f, g, h);
    // D + D: 2 f, 2 g, 2 h;  D * D: f^2, 2 f g, 2 (g g^T + f h)
    check_close(F_sum(q), 2. * f);
    check_close(F_prod(q), f * f);
    for (int i = 0; i < n; ++i) {
      CAPTURE(i);
      const auto ui = static_cast<std::size_t>(i);
      check_close(g_sum(n * q + i), 2. * g[ui]);
      check_close(g_prod(n * q + i), 2. * f * g[ui]);
      for (int j = 0; j < n; ++j) {
        CAPTURE(j);
        const auto uij = static_cast<std::size_t>(i * n + j);
        check_close(h_sum(n * n * q + i * n + j), 2. * h[uij]);
        check_close(h_prod(n * n * q + i * n + j),
                    2. * (g[ui] * g[static_cast<std::size_t>(j)] + f * h[uij]));
      }
    }
  }
}

TEST_CASE("eval_variational_at_points on the debug device: gradients (EnergySemiImplicit)",
          "[Coefficients][variational][device]") {
  use_debug_device();
  using EVAL = EnergySemiImplicitVariationalEval;
  constexpr int d = EVAL::dimension;
  constexpr double time = 3.;
  mfem::Vector values(npts), explicit_values(npts), auxiliary(npts);
  mfem::Vector gradients(d * npts), explicit_gradients(d * npts), auxiliary_gradients(d * npts);
  mfem::Vector coordinates(d * npts);
  for (int q = 0; q < npts; ++q) {
    values(q) = 0.1 + 0.2 * q;
    explicit_values(q) = 0.3 * q;
    auxiliary(q) = 1. - 0.1 * q;
    for (int k = 0; k < d; ++k) {
      gradients(d * q + k) = 0.5 * q - k;
      explicit_gradients(d * q + k) = 0.2 + 0.1 * q * (k + 1);
      auxiliary_gradients(d * q + k) = -0.3 + 0.25 * q + 0.5 * k;
      coordinates(d * q + k) = 1. + 0.5 * q + k;
    }
  }
  mfem::Vector F(npts), df_dvar(npts), df_dgrad(d * npts), d2f_dvar2(npts);
  mfem::Vector d2f_dvar_dgrad(d * npts), d2f_dgrad2(d * d * npts);
  VariationalBatch outputs;
  outputs.F = &F;
  outputs.df_dvar = &df_dvar;
  outputs.df_dgrad = &df_dgrad;
  outputs.d2f_dvar2 = &d2f_dvar2;
  outputs.d2f_dvar_dgrad = &d2f_dvar_dgrad;
  outputs.d2f_dgrad2 = &d2f_dgrad2;
  eval_variational_at_points(EVAL{}, npts, time, values, explicit_values, auxiliary, gradients,
                             explicit_gradients, auxiliary_gradients, coordinates, outputs);

  // Results written on the device: copied back to the host before being read
  for (mfem::Vector* v : {&F, &df_dvar, &df_dgrad, &d2f_dvar2, &d2f_dvar_dgrad, &d2f_dgrad2,
                          &values, &explicit_values, &auxiliary, &gradients, &explicit_gradients,
                          &auxiliary_gradients, &coordinates}) {
    v->HostRead();
  }
  for (int q = 0; q < npts; ++q) {
    CAPTURE(q);
    // Same structure evaluated on the host at point q
    double f = 0.;
    std::vector<double> f_u(1), f_p(d), f_uu(1), f_up(d), f_pp(d * d);
    VariationalOutputs out;
    out.F = &f;
    out.df_dvar = f_u.data();
    out.df_dgrad = f_p.data();
    out.d2f_dvar2 = f_uu.data();
    out.d2f_dvar_dgrad = f_up.data();
    out.d2f_dgrad2 = f_pp.data();
    const double u = values(q), ue = explicit_values(q), c = auxiliary(q);
    const std::vector<double> g{gradients(d * q), gradients(d * q + 1)};
    const std::vector<double> ge{explicit_gradients(d * q), explicit_gradients(d * q + 1)};
    const std::vector<double> gc{auxiliary_gradients(d * q), auxiliary_gradients(d * q + 1)};
    const std::vector<double> x{coordinates(d * q), coordinates(d * q + 1)};
    EVAL{}(&u, &ue, &c, g.data(), ge.data(), gc.data(), x.data(), time, out);

    check_close(F(q), f);
    check_close(df_dvar(q), f_u[0]);
    check_close(d2f_dvar2(q), f_uu[0]);
    for (int k = 0; k < d; ++k) {
      CAPTURE(k);
      check_close(df_dgrad(d * q + k), f_p[static_cast<std::size_t>(k)]);
      check_close(d2f_dvar_dgrad(d * q + k), f_up[static_cast<std::size_t>(k)]);
      for (int l = 0; l < d; ++l) {
        CAPTURE(l);
        check_close(d2f_dgrad2(d * d * q + k * d + l), f_pp[static_cast<std::size_t>(k * d + l)]);
      }
    }
    // Explicit term of the energy: 1e-3 T X u, T = time, X first coordinate
    check_close(df_dvar(q), 1.5 * (gc[0] * gc[0] + gc[1] * gc[1]) + 1.e-3 * time * x[0]);
  }
}
